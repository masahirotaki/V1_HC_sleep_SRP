function results = run_assembly_reactivation(subject, date, stimulus_name, event_rec, cfg, opts)
% RUN_ASSEMBLY_REACTIVATION  Reactivation of awake stimulus-element cell assemblies during
% ripples (sleep, or awake grey-screen ripples), and which element (A/B/C/D) each 20 ms bin
% carries.
%
%   results = run_assembly_reactivation('M00069', '20250929', 'GAVNIK_ABCD_1', 'Sleep_Box_2', project_config())
%   results = run_assembly_reactivation('M00069', '20250929', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_1', ...)  % awake ripples
%
% Needs run_sleep_detection on sleep recordings and run_abcd_awake_replay on the stimulus
% recording (its awake ripples are left out of the awake model and are the awake events).
%
% 1. Awake model (all sequences of STIMULUS_NAME, bins overlapping awake ripples +-50 ms left out):
%    spike counts in opts.bin (20 ms) bins from -0.6 s to sequence + 0.2 s around A onset;
%    each neuron z-scored over time, then the common mode (mean z across neurons in each bin)
%    subtracted; assemblies = PCA -> Marchenko-Pastur -> ICA (assembly_templates). Kept:
%    stimulus-element assemblies = activation during A-D bins > grey bins across sequences
%    (sign-rank, FDR opts.fdr) AND preferred class A-D. Element readouts: shrinkage LDA on the
%    kept assemblies' activations and a count-conditioned multinomial on the raw counts
%    (5 classes A/B/C/D/grey, labelled bins).
% 2. Pass 1 (per ripple; window [onset, max(offset, onset + 100 ms)] cut into 20 ms bins):
%    total activation T(t) = sum of the kept assemblies' activation strength; statistic =
%    max over bins (and, secondary, sum over bins). Null: neuron identities permuted with one
%    permutation for the whole ripple (opts.n_shuffle); significant if p < opts.alpha. The
%    shuffle picks its own maximum bin, so the number of bins is accounted for.
%    Z-scoring of event bins uses the event's own epoch: SWS / quiet-wake 20 ms bins for
%    sleep ripples (by state), grey-screen bins for awake ripples; common mode removed.
%    Also the ripple-triggered average of T(t) and of each assembly's activation (+-0.5 s,
%    aligned to ripple onset by default, opts.align = 'peak' for the peak), real vs shuffled.
%    Activation can be negative: R = sum over pairs i ~= j of w_i w_j z_i z_j, so it is < 0
%    when neurons with opposite-sign weights are co-active (or same-sign ones anti-active).
% 3. Pass 2 (significant ripples): element per 20 ms bin by LDA and by multinomial; null:
%    the same ripples after cell-ID shuffles (opts.n_shuffle_pass2 per ripple). The classes
%    include grey (the awake grey-screen pattern = no stimulus assembly active), and the
%    pass-1 statistic is the max over bins, so most bins of a significant ripple can be grey;
%    'active bins' (total activation > that ripple's shuffle 95% of the max) are also reported.
%    Order score: Spearman rho of element index vs bin index over the non-grey bins.
%    Examples: the opts.n_examples significant ripples with the largest max / threshold ratio,
%    per-assembly and total activation from -500 to +500 ms around ripple onset
%    (opts.example_window), with the decoded element per bin.
% opts.common_mode = false skips the common-mode removal everywhere (model and events);
% opts.suffix (e.g. '_nocm') keeps such variant outputs apart from the default ones.
%
% Outputs: <event_rec>\replay\replay_<event_rec>_assembly_reactivation.mat / .fig / .pdf
%          and replay_<event_rec>_assembly_reactivation_examples.fig / .pdf
%          (awake model also in <stimulus>\Decoding\decoding_<stimulus>_assembly_model.mat)

if nargin < 5 || isempty(cfg); cfg = project_config(); end
if nargin < 6; opts = struct(); end
def.populations = struct('name', {'CA1_all', 'V1_good'}, 'region', {"CA1", "V1"}, 'good_only', {false, true});
def.bin = 0.02;
def.grey = [-0.6 -0.45];
def.min_spikes = 20;
def.fdr = 0.05;
def.n_shuffle = 1000;
def.n_shuffle_pass2 = 20;
def.alpha = 0.05;
def.states = {'SWS', 'quietWake'};
def.peri = [-0.5 0.5];
def.n_shuffle_peri = 20;
def.grey_margin = 0.3;
def.align = 'onset';          % ripple-triggered average aligned to ripple 'onset' or 'peak'
def.n_examples = 8;
def.example_window = [-0.5 0.5];   % s around ripple onset shown in the example figure
def.common_mode = true;      % subtract the across-neuron mean z per bin (false = per-neuron z only)
def.suffix = '';             % appended to output names (e.g. '_nocm' for a variant run)
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
awake_events = strcmp(event_rec, stimulus_name);
tag = sprintf('%s %s %s (assemblies from %s)', subject, date, strrep(event_rec, '_', ' '), strrep(stimulus_name, '_', ' '));

%% awake model
spaths = get_session_paths(subject, date, stimulus_name, cfg);
sunits = load_session_clusters(spaths, get_probe_layout(spaths, cfg), cfg);
trials = get_abcd_trials(spaths);
AR = load(fullfile(spaths.analysis_dir, 'replay', sprintf('replay_%s_results.mat', stimulus_name)));
awake_rip = AR.results.ripples;
model = build_model(sunits, trials, awake_rip, opts);

%% events
if awake_events
    units = sunits;
    rip = awake_rip;
    C = classify_ripple_clusters(rip.onset(:), rip.offset(:));
    pk = rip.peaktimes(C.order); pk = pk(:);
    ctr = mean(C.window, 2);
    grey_ok = false(size(ctr));
    for k = 1:numel(trials.seq_onset) - 1
        grey_ok = grey_ok | (ctr >= trials.seq_offset(k) + opts.grey_margin & C.window(:, 2) < trials.seq_onset(k + 1));
    end
    state = repmat("awake_grey", numel(ctr), 1);
    state(~grey_ok) = "stimulus";
    epochs.awake_grey = grey_intervals(trials, opts.grey_margin);
    use_states = "awake_grey";
    out_dir = fullfile(spaths.analysis_dir, 'replay');
else
    paths = get_session_paths(subject, date, event_rec, cfg);
    units = load_session_clusters(paths, get_probe_layout(paths, cfg), cfg);
    R = load(fullfile(paths.detection_dir, 'extracted_ripple_events.mat')); rip = R.ripples;
    B = load(fullfile(paths.output_dir, 'behavioural_state.mat')); states = B.behavioural_state;
    C = classify_ripple_clusters(rip.onset(:), rip.offset(:));
    pk = rip.peaktimes(C.order); pk = pk(:);
    state = repmat("other", numel(C.onset), 1);
    for s = string(opts.states)
        iv = states.(s);
        for k = 1:size(iv, 1); state(C.onset >= iv(k, 1) & C.onset < iv(k, 2)) = s; end
        epochs.(s) = iv;
    end
    use_states = string(opts.states);
    out_dir = fullfile(paths.analysis_dir, 'replay');
end
check_output_path(out_dir, cfg);
scored = ismember(state, use_states);
fprintf('%s: %d ripples, %d scored (%s)\n', event_rec, numel(C.onset), sum(scored), strjoin(use_states, '/'));

%% passes per population
for p = 1:numel(model)
    M = model(p);
    out = struct('name', M.name, 'n_assemblies', size(M.W_all, 2), 'n_kept', size(M.W, 2), 'preferred', M.pref, ...
        'kept', M.keep, 'p_stim', M.p_stim);
    if isempty(M.W)
        fprintf('  %s: no stimulus-element assemblies (of %d) - skipped\n', M.name, size(M.W_all, 2));
        pops(p) = fill_empty(out); %#ok<AGROW>
        continue
    end
    present = ismember(M.cluster_id, units.cluster_id);
    if ~all(present); warning('%d %s template units are missing in %s', sum(~present), M.name, event_rec); end
    % epoch statistics (per state) for z-scoring
    for s = use_states
        [stat.(s).mu, stat.(s).sd] = epoch_stats(units, M.cluster_id, epochs.(s), opts.bin, C.window);
    end
    ev = find(scored);
    nev = numel(ev);
    [obs_max, obs_sum, p_max, p_sum, nbins, null95] = deal(nan(nev, 1));
    Zc = cell(nev, 1); Cc = cell(nev, 1);
    w = collect_spike_windows(units.spike_times, units.spike_id, M.cluster_id, C.window(ev, 1), ...
        C.window(ev, 1) + ceil((C.window(ev, 2) - C.window(ev, 1)) / opts.bin - 1e-9) * opts.bin);
    for k = 1:nev
        nb = ceil((C.window(ev(k), 2) - C.window(ev(k), 1)) / opts.bin - 1e-9);
        cnt = zeros(nb, numel(M.cluster_id));
        if ~isempty(w{k})
            b = min(floor(w{k}(:, 1) / opts.bin) + 1, nb);
            cnt = accumarray([b w{k}(:, 2)], 1, size(cnt));
        end
        s = state(ev(k));
        Z = (cnt - stat.(s).mu) ./ stat.(s).sd;
        Zc{k} = Z - opts.common_mode * mean(Z, 2);
        Cc{k} = cnt;
        nbins(k) = nb;
    end
    W = M.W; nsh = opts.n_shuffle; seed0 = opts.seed + 100 * p; opts_alpha = opts.alpha;
    parfor k = 1:nev
        Z = Zc{k};
        T = total_activation(Z, W);
        nm = zeros(nsh, 1); ns = nm;
        rs = RandStream('mt19937ar', 'Seed', seed0 + k);
        for q = 1:nsh
            Tq = total_activation(Z(:, randperm(rs, size(Z, 2))), W);
            nm(q) = max(Tq); ns(q) = sum(Tq);
        end
        obs_max(k) = max(T); obs_sum(k) = sum(T);
        null95(k) = prctile(nm, 100 * (1 - opts_alpha));        % shuffle threshold for the max
        p_max(k) = (1 + sum(nm >= obs_max(k))) / (1 + nsh);
        p_sum(k) = (1 + sum(ns >= obs_sum(k))) / (1 + nsh);
    end
    sig = p_max < opts.alpha;
    % pass 2: element per bin in significant ripples, real and cell-ID shuffled
    % (5 classes A/B/C/D/grey: grey = the awake grey-screen pattern, i.e. no stimulus assembly
    % active; 'active' bins = bins whose total activation exceeds the ripple's shuffle threshold)
    [lda_cls, mn_cls] = deal(cell(nev, 1));
    sh_lda = []; sh_mn = []; act_lda = []; act_mn = [];
    rng(opts.seed + 7 * p);
    for k = find(sig)'
        Ra = assembly_activation(Zc{k}, W);
        lda_cls{k} = apply_decoder(M.lda, Ra);
        mn_cls{k} = apply_decoder(M.mnb, Cc{k});
        a = sum(Ra, 2) > null95(k);
        act_lda = [act_lda; lda_cls{k}(a)]; act_mn = [act_mn; mn_cls{k}(a)]; %#ok<AGROW>
        for q = 1:opts.n_shuffle_pass2
            pm = randperm(size(Zc{k}, 2));
            sh_lda = [sh_lda; apply_decoder(M.lda, assembly_activation(Zc{k}(:, pm), W))]; %#ok<AGROW>
            sh_mn = [sh_mn; apply_decoder(M.mnb, Cc{k}(:, pm))]; %#ok<AGROW>
        end
    end
    real_lda = cell2mat(lda_cls(sig)); real_mn = cell2mat(mn_cls(sig));
    frac = @(x) accumarray(x(:), 1, [5 1])' / max(numel(x), 1);
    out.ripple = table(ev, C.onset(ev), C.window(ev, 2), C.type(ev), state(ev), nbins, obs_max, null95, p_max, obs_sum, p_sum, sig, ...
        'VariableNames', {'ripple_index', 'onset', 'window_end', 'type', 'state', 'n_bins', 'max_T', 'null95_max', 'p_max', ...
        'sum_T', 'p_sum', 'significant'});
    out.element_lda = lda_cls; out.element_multinomial = mn_cls;
    out.class_frac = struct('lda', frac(real_lda), 'lda_shuffle', frac(sh_lda), ...
        'multinomial', frac(real_mn), 'multinomial_shuffle', frac(sh_mn), ...
        'lda_active', frac(act_lda), 'multinomial_active', frac(act_mn), 'n_active_bins', numel(act_lda));
    out.order_rho = cellfun(@order_rho, mn_cls(sig));
    % event-triggered average of total and per-assembly activation (main state), aligned to
    % ripple onset (opts.align = 'onset') or ripple peak ('peak')
    s0 = use_states(1);
    if strcmp(opts.align, 'peak'); t_al = pk(state == s0); else; t_al = C.onset(state == s0); end
    [out.peri_t, out.peri_T, out.peri_T_shuffle, out.peri_R, out.peri_R_shuffle] = peri_average(units, M, stat.(s0), t_al, W, opts);
    out.align = opts.align;
    % example significant ripples: activation time courses and decoded element per bin
    out.examples = example_events(units, M, stat, state, ev, sig, obs_max ./ null95, null95, C, opts);
    out.alpha = opts.alpha;
    fprintf(['  %s: %d/%d assemblies kept (prefer A-D: %s); significant ripples %d/%d (%.1f%%, binomial p = %.2g); ' ...
        'by sum %.1f%%\n'], M.name, size(W, 2), size(M.W_all, 2), mat2str(accumarray(M.pref(M.keep), 1, [4 1])'), sum(sig), nev, ...
        100 * mean(sig), 1 - binocdf(sum(sig) - 1, nev, opts.alpha), 100 * mean(p_sum < opts.alpha));
    fprintf('     pass 2 multinomial classes A B C D grey: real %s | shuffle %s\n', ...
        mat2str(round(100 * out.class_frac.multinomial)), mat2str(round(100 * out.class_frac.multinomial_shuffle)));
    pops(p) = out; %#ok<AGROW>
end

results = struct('subject', subject, 'date', date, 'stimulus', stimulus_name, 'event_recording', event_rec, ...
    'opts', opts, 'clusters', C, 'state', state, 'pop', pops);
model_out = rmfield(model, 'Zawake');
fig = [plot_reactivation(results, tag), plot_examples(results, tag)];
if ~opts.save; return; end
safe_save(fullfile(spaths.analysis_dir, 'Decoding', sprintf('decoding_%s_assembly_model%s.mat', stimulus_name, opts.suffix)), ...
    struct('model', model_out), cfg);
base = sprintf('replay_%s_assembly_reactivation%s', event_rec, opts.suffix);
safe_save(fullfile(out_dir, [base '.mat']), struct('results', results), cfg);
save_vector_figures(fig, {base, [base '_examples']}, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function model = build_model(units, trials, awake_rip, opts)
seq = round(trials.seq_duration, 2);
edges = opts.grey(1):opts.bin:seq + 0.2 + 1e-9;
tc = edges(1:end - 1) + opts.bin / 2;
nb = numel(tc);
n_el = numel(trials.element_names);
lab = zeros(1, nb);
for e = 1:n_el
    lab(tc >= trials.element_start(e) & tc < trials.element_start(e) + trials.element_duration) = e;
end
lab(tc >= opts.grey(1) & tc < opts.grey(2)) = n_el + 1;
n_seq = numel(trials.seq_onset);
abs_t = trials.seq_onset + tc;                                  % seq x bin
bad = false(size(abs_t));
for k = 1:numel(awake_rip.onset)
    bad = bad | (abs_t + opts.bin / 2 > awake_rip.onset(k) - 0.05 & abs_t - opts.bin / 2 < awake_rip.offset(k) + 0.05);
end
for p = 1:numel(opts.populations)
    P = opts.populations(p);
    ids = units.cluster_id(units.region == P.region & (units.is_good | ~P.good_only));
    X = bin_trial_spikes(units.spike_times, units.spike_id, ids, trials.seq_onset, edges);   % seq x neuron x bin
    keepn = squeeze(sum(X, [1 3])) >= opts.min_spikes;
    X = X(:, keepn, :); ids = ids(keepn);
    Xb = reshape(permute(X, [3 1 2]), nb * n_seq, []);
    good = ~reshape(bad', [], 1);
    lab_b = repmat(lab', n_seq, 1);
    seq_b = repelem((1:n_seq)', nb);
    mu = mean(Xb(good, :), 1); sd = std(Xb(good, :), 0, 1); sd(sd == 0) = 1;
    Z = (Xb - mu) ./ sd; Z = Z - opts.common_mode * mean(Z, 2);
    A = assembly_templates(Z(good, :), 'ica', opts.seed);
    R = assembly_activation(Z, A.W);
    k = size(A.W, 2);
    prof = zeros(k, n_el + 1); p_stim = ones(k, 1);
    for c = 1:n_el + 1; prof(:, c) = mean(R(good & lab_b == c, :), 1)'; end
    [~, pref] = max(prof, [], 2);
    for j = 1:k
        st = accumarray(seq_b(good & lab_b >= 1 & lab_b <= n_el), R(good & lab_b >= 1 & lab_b <= n_el, j), [n_seq 1], @mean, NaN);
        gr = accumarray(seq_b(good & lab_b == n_el + 1), R(good & lab_b == n_el + 1, j), [n_seq 1], @mean, NaN);
        ok = ~isnan(st) & ~isnan(gr);
        if sum(ok) > 10; p_stim(j) = signrank(st(ok), gr(ok), 'tail', 'right'); end
    end
    sig = fdr(p_stim, opts.fdr);
    keep = sig & pref <= n_el;
    L = good & lab_b > 0;
    model(p).name = P.name; %#ok<AGROW>
    model(p).cluster_id = ids;
    model(p).mu = mu; model(p).sd = sd;
    model(p).W_all = A.W; model(p).W = A.W(:, keep);
    model(p).profile = prof; model(p).pref = pref; model(p).p_stim = p_stim; model(p).keep = keep;
    model(p).bin = opts.bin;
    if any(keep)
        model(p).lda = fit_decoder(R(L, keep), lab_b(L), n_el + 1, 'lda');
    else
        model(p).lda = [];
    end
    model(p).mnb = fit_decoder(Xb(L, :), lab_b(L), n_el + 1, 'multinomial');
    model(p).Zawake = [];
end
end

function T = total_activation(Z, W)
T = sum((Z * W).^2 - (Z.^2) * (W.^2), 2);
end

function [mu, sd] = epoch_stats(units, ids, iv, bin, excl)
% per-neuron mean / SD of 20 ms counts inside the intervals iv, ripple windows excluded
in = ismember(units.spike_id, ids);
st = units.spike_times(in);
[~, col] = ismember(units.spike_id(in), ids);
s1 = zeros(1, numel(ids)); s2 = s1; nb = 0;
for k = 1:size(iv, 1)
    e = iv(k, 1):bin:iv(k, 2);
    if numel(e) < 2; continue; end
    m = st >= e(1) & st < e(end);
    if any(m)
        c = accumarray([discretize(st(m), e), col(m)], 1, [numel(e) - 1, numel(ids)]);
    else
        c = zeros(numel(e) - 1, numel(ids));
    end
    ctr = e(1:end - 1) + bin / 2;
    drop = any(ctr' >= excl(:, 1)' - 0.05 & ctr' <= excl(:, 2)' + 0.05, 2);
    c = c(~drop, :);
    s1 = s1 + sum(c, 1); s2 = s2 + sum(c.^2, 1); nb = nb + size(c, 1);
end
mu = s1 / max(nb, 1);
sd = sqrt(max(s2 / max(nb, 1) - mu.^2, 0) * nb / max(nb - 1, 1));
sd(sd == 0) = 1;
end

function iv = grey_intervals(trials, margin)
iv = [trials.seq_offset(1:end - 1) + margin, trials.seq_onset(2:end)];
iv = iv(iv(:, 2) > iv(:, 1), :);
end

function [t, Tm, Ts, Rm, Rs] = peri_average(units, M, stat, t_al, W, opts)
% mean total (T) and per-assembly (R) activation around the alignment times, real and with
% the cell-ID shuffle (one permutation per event, opts.n_shuffle_peri of them)
edges = opts.peri(1):opts.bin:opts.peri(2) + 1e-9;
t = edges(1:end - 1) + opts.bin / 2;
nb = numel(t); k_as = size(W, 2);
Tm = nan(1, nb); Ts = nan(1, nb); Rm = nan(nb, k_as); Rs = nan(nb, k_as);
if isempty(t_al); return; end
w = collect_spike_windows(units.spike_times, units.spike_id, M.cluster_id, t_al + edges(1), t_al + edges(end));
acc = zeros(nb, k_as); accs = zeros(nb, k_as);
rng(opts.seed);
for k = 1:numel(t_al)
    cnt = zeros(nb, numel(M.cluster_id));
    if ~isempty(w{k}); cnt = accumarray([min(floor(w{k}(:, 1) / opts.bin) + 1, nb), w{k}(:, 2)], 1, size(cnt)); end
    Z = (cnt - stat.mu) ./ stat.sd; Z = Z - opts.common_mode * mean(Z, 2);
    acc = acc + assembly_activation(Z, W);
    for q = 1:opts.n_shuffle_peri
        accs = accs + assembly_activation(Z(:, randperm(size(Z, 2))), W) / opts.n_shuffle_peri;
    end
end
Rm = acc / numel(t_al); Rs = accs / numel(t_al);
Tm = sum(Rm, 2)'; Ts = sum(Rs, 2)';
end

function X = example_events(units, M, stat, state, ev, sig, score, null95, C, opts)
% the opts.n_examples significant ripples with the largest max_T / shuffle threshold:
% per-assembly activation and decoded element per 20 ms bin in opts.example_window around
% the ripple onset (bins aligned to the onset; the window is extended if the ripple is longer)
X = struct('ripple_index', {}, 'onset', {}, 'window', {}, 'state', {}, 'type', {}, 't', {}, 'R', {}, 'T', {}, ...
    'null95_max', {}, 'lda', {}, 'multinomial', {}, 'other_ripples', {});
i = find(sig);
[~, o] = sort(score(i), 'descend');
i = i(o(1:min(opts.n_examples, numel(i))));
for j = 1:numel(i)
    k = i(j); r = ev(k);
    nb_pre = round(-opts.example_window(1) / opts.bin);
    nb_win = ceil((C.window(r, 2) - C.window(r, 1)) / opts.bin - 1e-9);
    nb = nb_pre + max(round(opts.example_window(2) / opts.bin), nb_win);
    t0 = C.window(r, 1) - nb_pre * opts.bin;
    w = collect_spike_windows(units.spike_times, units.spike_id, M.cluster_id, t0, t0 + nb * opts.bin); w = w{1};
    cnt = zeros(nb, numel(M.cluster_id));
    if ~isempty(w); cnt = accumarray([min(floor(w(:, 1) / opts.bin) + 1, nb), w(:, 2)], 1, size(cnt)); end
    s = state(r);
    Z = (cnt - stat.(s).mu) ./ stat.(s).sd; Z = Z - opts.common_mode * mean(Z, 2);
    R = assembly_activation(Z, M.W);
    X(j).ripple_index = r; X(j).onset = C.window(r, 1); X(j).window = C.window(r, :) - C.window(r, 1);
    X(j).state = s; X(j).type = C.type(r);
    X(j).t = ((0:nb - 1) + 0.5) * opts.bin - nb_pre * opts.bin;     % bin centres re onset
    X(j).R = R; X(j).T = sum(R, 2); X(j).null95_max = null95(k);
    X(j).lda = apply_decoder(M.lda, R); X(j).multinomial = apply_decoder(M.mnb, cnt);
    near = find(C.window(:, 2) > t0 & C.window(:, 1) < t0 + nb * opts.bin & (1:size(C.window, 1))' ~= r);
    X(j).other_ripples = C.window(near, :) - C.window(r, 1);
end
end

function r = order_rho(x)
x = x(:); i = find(x <= 4);
if numel(i) < 2 || numel(unique(x(i))) < 2; r = NaN; return; end
r = corr(i, x(i), 'type', 'Spearman');
end

function s = fdr(p, q)
[ps, o] = sort(p(:)); m = numel(ps);
k = find(ps <= (1:m)' / m * q, 1, 'last');
s = false(size(p)); if ~isempty(k); s(o(1:k)) = true; end
end

function out = fill_empty(out)
out.ripple = table(); out.element_lda = {}; out.element_multinomial = {};
out.class_frac = struct('lda', nan(1, 5), 'lda_shuffle', nan(1, 5), 'multinomial', nan(1, 5), 'multinomial_shuffle', nan(1, 5), ...
    'lda_active', nan(1, 5), 'multinomial_active', nan(1, 5), 'n_active_bins', 0);
out.order_rho = []; out.peri_t = []; out.peri_T = []; out.peri_T_shuffle = []; out.peri_R = []; out.peri_R_shuffle = [];
out.align = ''; out.examples = []; out.alpha = NaN;
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function c = element_colours()
c = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655; 0.8 0.8 0.8];
end

function f = plot_reactivation(r, tag)
np = numel(r.pop);
blue = [0.165 0.471 0.839]; grey = [0.55 0.55 0.55]; orange = [0.922 0.408 0.204];
ccol = element_colours();
cls = {'A', 'B', 'C', 'D', 'grey'};
f = figure('Color', 'w', 'Position', [30 30 2100 380 * np]);
tl = tiledlayout(f, np, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - stimulus-element assembly reactivation in ripples (20 ms bins; pass 1: max over bins of total ' ...
    'activation vs cell-ID shuffle; pass 2: element per bin, classes A/B/C/D/grey)'], tag), 'FontSize', 10.5);
for p = 1:np
    P = r.pop(p);
    pname = sprintf('%s: %d of %d assemblies kept', strrep(P.name, '_', ' '), P.n_kept, P.n_assemblies);
    if isempty(P.ripple) || height(P.ripple) == 0
        ax = nexttile(tl, [1 5]); axis(ax, 'off');
        text(ax, 0.5, 0.5, [pname ' - no stimulus-element assemblies, nothing scored'], 'HorizontalAlignment', 'center');
        continue
    end
    T = P.ripple;
    al = P.align; if isempty(al); al = 'peak'; end

    % 1. % significant
    ax = nexttile(tl); hold(ax, 'on');
    groups = {'all', 'solo', 'cluster', 'intermediate'};
    v = zeros(1, 4); n = v;
    for g = 1:4
        m = true(height(T), 1); if g > 1; m = T.type == groups{g}; end
        v(g) = 100 * mean(T.significant(m)); n(g) = sum(m);
    end
    b = bar(ax, v, 'FaceColor', blue, 'EdgeColor', 'none');
    text(ax, b.XEndPoints, v + max(v) * 0.03, compose('%.0f%%', v), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 7);
    yline(ax, 100 * P.alpha, ':', 'chance', 'FontSize', 7);
    ylim(ax, [0 max([v 100 * P.alpha]) * 1.2]);
    set(ax, 'XTick', 1:4, 'XTickLabel', strcat(groups, " (", string(n), ")")); xtickangle(ax, 25);
    ylabel(ax, '% ripples with significant reactivation');
    title(ax, sprintf('%s\nbinomial p = %.2g', pname, 1 - binocdf(sum(T.significant) - 1, height(T), P.alpha)), ...
        'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);

    % 2. event-triggered total activation
    ax = nexttile(tl); hold(ax, 'on');
    if ~isempty(P.peri_t)
        plot(ax, 1000 * P.peri_t, P.peri_T_shuffle, '-', 'Color', grey, 'LineWidth', 1.5, 'DisplayName', 'cell-ID shuffled');
        plot(ax, 1000 * P.peri_t, P.peri_T, '-', 'Color', blue, 'LineWidth', 2, 'DisplayName', 'real');
        xline(ax, 0, 'k:', 'HandleVisibility', 'off'); yline(ax, 0, 'k-', 'HandleVisibility', 'off');
        legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
    end
    xlabel(ax, sprintf('time from ripple %s (ms)', al)); ylabel(ax, 'total stimulus-assembly activation');
    title(ax, sprintf('ripple-%s-triggered average, total (main state)', al), 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);

    % 3. event-triggered per-assembly activation, coloured by preferred element
    ax = nexttile(tl); hold(ax, 'on');
    if isfield(P, 'peri_R') && ~isempty(P.peri_R)
        pref = P.preferred(P.kept);
        for j = 1:size(P.peri_R, 2)
            plot(ax, 1000 * P.peri_t, P.peri_R(:, j), '-', 'Color', ccol(pref(j), :), 'LineWidth', 1.2);
        end
        plot(ax, 1000 * P.peri_t, mean(P.peri_R_shuffle, 2), '-', 'Color', grey, 'LineWidth', 1.5);
        xline(ax, 0, 'k:'); yline(ax, 0, 'k-');
        for e = 1:4
            plot(ax, NaN, NaN, '-', 'Color', ccol(e, :), 'LineWidth', 1.5, 'DisplayName', sprintf('prefers %s', cls{e}));
        end
        hh = findobj(ax, 'Type', 'line', '-regexp', 'DisplayName', 'prefers');
        legend(ax, flipud(hh), 'Location', 'best', 'Box', 'off', 'FontSize', 7);
    end
    xlabel(ax, sprintf('time from ripple %s (ms)', al)); ylabel(ax, 'activation per assembly');
    title(ax, 'per assembly (grey = shuffled mean)', 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);

    % 4. pass 2 composition: all bins of significant ripples vs active bins only
    ax = nexttile(tl); hold(ax, 'on');
    cf = P.class_frac;
    if ~isfield(cf, 'lda_active'); cf.lda_active = nan(1, 5); cf.multinomial_active = nan(1, 5); cf.n_active_bins = 0; end
    F = [cf.lda; cf.lda_shuffle; cf.lda_active; cf.multinomial; cf.multinomial_shuffle; cf.multinomial_active]' * 100;
    bb = bar(ax, F, 'BarWidth', 0.9, 'EdgeColor', 'none');
    fc = [blue; blue * 0.35 + 0.65; blue * 0.6; orange; orange * 0.35 + 0.65; orange * 0.6];
    for j = 1:6; bb(j).FaceColor = fc(j, :); end
    set(ax, 'XTick', 1:5, 'XTickLabel', cls);
    ylabel(ax, '% of bins');
    legend(ax, {'LDA, all bins', 'LDA, shuffled', 'LDA, active bins', 'multinomial, all bins', 'multinomial, shuffled', ...
        'multinomial, active bins'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 6.5);
    title(ax, sprintf('pass 2: %d sig. ripples; %d active bins', sum(T.significant), ...
        cf.n_active_bins), 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);

    % 5. element per bin, all significant ripples
    ax = nexttile(tl);
    E = P.element_multinomial(T.significant);
    mx = max([cellfun(@numel, E); 1]);
    I = nan(numel(E), mx);
    for k = 1:numel(E); I(k, 1:numel(E{k})) = E{k}; end
    if ~isempty(I)
        imagesc(ax, I, 'AlphaData', ~isnan(I), [0.5 5.5]); colormap(ax, ccol);
        colorbar(ax, 'Ticks', 1:5, 'TickLabels', cls);
    end
    xlabel(ax, '20 ms bin within ripple'); ylabel(ax, 'significant ripple');
    rr = P.order_rho(~isnan(P.order_rho));
    title(ax, sprintf('multinomial per bin; order rho %.2f (n = %d)', mean(rr), numel(rr)), ...
        'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
end
end

function f = plot_examples(r, tag)
% one figure: for each population, up to 8 example significant ripples; per event the
% per-assembly activation (thin, coloured by preferred element), total activation (black),
% the ripple's shuffle 95% threshold for the max (dashed), ripple windows (grey), and below
% the element decoded in each bin (LDA on assemblies, multinomial on counts)
ccol = element_colours(); cls = {'A', 'B', 'C', 'D', 'grey'};
pops = r.pop(arrayfun(@(P) ~isempty(P.examples), r.pop));
f = figure('Color', 'w', 'Position', [20 20 2000 520 * max(numel(pops), 1)]);
ncol = 4;
nrow_ev = sum(arrayfun(@(P) ceil(numel(P.examples) / ncol), pops));
tl = tiledlayout(f, max(3 * nrow_ev, 1), ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - example significant ripples: activation per assembly (colour = preferred element), total ' ...
    '(black), shuffle 95%% of the max (dashed), ripples (grey); strips = decoded element per 20 ms bin: A blue, B orange, C green, D purple, grey'], tag), 'FontSize', 10);
if isempty(pops)
    ax = nexttile(tl); axis(ax, 'off'); text(ax, 0.5, 0.5, 'no significant ripples', 'HorizontalAlignment', 'center'); return
end
row0 = 0;
for p = 1:numel(pops)
    P = pops(p); X = P.examples; pref = P.preferred(P.kept);
    for k = 1:numel(X)
        rr = row0 + floor((k - 1) / ncol); cc = mod(k - 1, ncol) + 1;
        x = X(k); tms = 1000 * x.t;
        ax = nexttile(tl, (3 * rr) * ncol + cc, [2 1]); hold(ax, 'on');
        yl = [min([x.R(:); x.T; 0]) max([x.R(:); x.T; x.null95_max])] .* [1.1 1.15];
        for q = [x.window; x.other_ripples]'
            patch(ax, 1000 * q([1 2 2 1])', yl([1 1 2 2]), [0.92 0.92 0.92], 'EdgeColor', 'none');
        end
        for j = 1:size(x.R, 2)
            stairs(ax, tms - 10, x.R(:, j), '-', 'Color', ccol(pref(j), :), 'LineWidth', 1);
        end
        stairs(ax, tms - 10, x.T, 'k-', 'LineWidth', 1.8);
        yline(ax, x.null95_max, 'k--'); yline(ax, 0, 'k:'); xline(ax, 0, 'k:');
        ylim(ax, yl); xlim(ax, [tms(1) - 10, tms(end) + 10]);
        set(ax, 'XTickLabel', []); style(ax);
        if cc == 1; ylabel(ax, 'activation'); end
        title(ax, sprintf('%s | ripple %d (%s, %s) | max T %.1f vs shuffle 95%% %.1f', strrep(P.name, '_', ' '), ...
            x.ripple_index, x.type, x.state, max(x.T(x.t >= 0 & x.t <= x.window(2))), x.null95_max), ...
            'FontWeight', 'normal', 'FontSize', 7.5);
        ax = nexttile(tl, (3 * rr + 2) * ncol + cc);
        imagesc(ax, tms, 1:2, [x.lda(:)'; x.multinomial(:)'], [0.5 5.5]); colormap(ax, ccol);
        hold(ax, 'on'); xline(ax, 0, 'k:'); xline(ax, 1000 * x.window(2), 'k:');
        set(ax, 'YTick', 1:2, 'YTickLabel', {'LDA', 'multin.'}, 'FontSize', 7, 'TickDir', 'out', 'Box', 'off');
        xlim(ax, [tms(1) - 10, tms(end) + 10]);
        xlabel(ax, 'ms from ripple onset');
    end
    row0 = row0 + ceil(numel(X) / ncol);
end
end
