function results = run_abcd_awake_replay(subject, date, stimulus_name, cfg, opts)
% RUN_ABCD_AWAKE_REPLAY  Awake ripples, MUA candidate events and sequence-likelihood
% replay (Huh et al. 2026) of the ABCD response order during a GAVNIK ABCD recording.
%
%   results = run_abcd_awake_replay('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Needs run_abcd_sequence_likelihood to have been run first (uses its folds and units).
%
% 1. Ripples: FindRipples_masa on the CA1 ripple channel (cfg.ripple: onset/offset z = 1,
%    peak z >= 5, 125-300 Hz).
% 2. Candidate events: detect_candidate_events_masa (VR_NPX_analysis; Davidson et al. 2009):
%    MUA of all CA1 units, 0.1 s Gaussian, z-scored over the recording; peak z >= 3,
%    edges at z = 0, >= 100 ms, merged if < 50 ms apart, >= 5 active neurons. Kept as
%    replay candidates if their ripple-band peak z >= opts.ripple_peak_min (5).
%    No movement exclusion: this head-fixed recording has no mobility trace.
% 3. Replay: CA1 spikes of each candidate scored with the CA1 ordered firing probability
%    matrix of the ABCD response (0-0.8 s after A; seq_order_matrix / seq_likelihood),
%    1000 cell-ID shuffles. Forward replay: Lseq above 95% of shuffles; reverse: below 5%
%    (chance 5% each, as in the paper). Candidates that overlap a training window are
%    scored with the matrix of the CV fold that left that sequence out (no circularity).
%    V1 spikes in the same events are scored with the V1 matrix as well.
%    Calibration: the cell-ID shuffle is liberal when spikes cluster in time (seen in the
%    awake grey windows), so each event is also compared with 20 random grey windows of
%    its own duration (no candidate events in them); "calibrated" replay must also be
%    beyond the 95th (5th) percentile of those windows' z.
% 4. Timing: ripples, candidates and significant replays relative to A onset (-0.5 to 1 s).
%
% Outputs (replay_* names) go to <analysis>/<subject>/<date>/<stimulus>/replay.

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
opts = fill_defaults(opts);
set(groot, 'defaultFigureVisible', cfg.plot.visible);

paths  = get_session_paths(subject, date, stimulus_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
trials = get_abcd_trials(paths);
out_dir = fullfile(paths.analysis_dir, 'replay');
check_output_path(out_dir, cfg);
tag = sprintf('%s %s %s', subject, date, strrep(stimulus_name, '_', ' '));
S = load(fullfile(paths.analysis_dir, 'Decoding', sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
seqlik = S.results;
fold = seqlik.fold;
win = seqlik.opts.window;

%% LFP and behaviour (no mobility in this recording -> no movement exclusion)
LFP = load_npx_lfp(paths.lf_bin_file, layout.chanMap(layout.ripple_channel), cfg);
ca1 = double(LFP.data(1, :));
tvec = LFP.tvec;
B = load(fullfile(paths.analysis_dir, 'extracted_behaviour.mat'));
fb = fieldnames(B);
B = B.(fb{1});
beh.sglxTime = unique(B.sglxTime(isfinite(B.sglxTime)))';
beh.tvec = beh.sglxTime;
beh.mobility_zscore = zeros(size(beh.sglxTime));
if isfield(B, 'mobility')
    warning('Behaviour has a mobility trace; this script does not use it for exclusion yet.');
end

%% 1. ripples
rip = FindRipples_masa(ca1(:), tvec(:), 'behaviour', beh, 'minDuration', cfg.ripple.minDuration, ...
    'durations', cfg.ripple.durations, 'frequency', LFP.fs, 'noise', [], 'passband', cfg.ripple.passband, ...
    'thresholds', cfg.ripple.thresholds, 'show', 'off', 'best_channel', layout.ripple_channel);
close(findall(groot, 'Type', 'figure'));
fprintf('Ripples: %d (%.1f / min)\n', numel(rip.peaktimes), numel(rip.peaktimes) / (range(tvec) / 60));

%% 2. MUA candidate events (all CA1 units)
ca1_all = units.cluster_id(units.region == "CA1");
in = ismember(units.spike_id, ca1_all);
mua_spikes = [double(units.spike_id(in)) units.spike_times(in)];
[cand, ~] = detect_candidate_events_masa(tvec, ca1, mua_spikes, beh, [], [], struct());
cand = struct_rows(cand);
cand.is_ripple = cand.ripple_peak(:) >= opts.ripple_peak_min;
cand.centre = (cand.onset(:) + cand.offset(:)) / 2;
fprintf('Candidate events (MUA z >= 3): %d, with ripple peak z >= %g: %d\n', numel(cand.onset), ...
    opts.ripple_peak_min, sum(cand.is_ripple));

%% 3. replay scoring (CA1 and V1 matrices, fold-aware)
ev = find(cand.is_ripple);
for p = 1:numel(opts.populations)
    name = opts.populations{p};
    P = seqlik.pop(strcmp({seqlik.pop.name}, name));
    ids = P.cluster_id;
    n = numel(ids);
    stim_w = collect_spike_windows(units.spike_times, units.spike_id, ids, trials.seq_onset + win(1), ...
        trials.seq_onset + win(2));
    M_all = seq_order_matrix(stim_w, n);
    M_fold = cell(max(fold), 1);
    for f = 1:max(fold)
        M_fold{f} = seq_order_matrix(stim_w(fold ~= f), n);
    end
    ev_w = collect_spike_windows(units.spike_times, units.spike_id, ids, cand.onset(ev), cand.offset(ev));
    % which matrix: the fold that left out any overlapping training sequence
    which_fold = zeros(numel(ev), 1);
    for k = 1:numel(ev)
        ov = find(cand.offset(ev(k)) > trials.seq_onset + win(1) & cand.onset(ev(k)) < trials.seq_onset + win(2), 1);
        if ~isempty(ov); which_fold(k) = fold(ov); end
    end
    R = repmat(struct('L', NaN, 'L_norm', NaN, 'z', NaN, 'pct', NaN, 'n_spikes', 0, 'n_neurons', 0), numel(ev), 1);
    ns = opts.n_shuffle;
    min_n = opts.min_neurons;
    parfor k = 1:numel(ev)
        w = ev_w{k}; %#ok<PFBNS>
        if numel(unique(w(:, 2))) < min_n; continue; end
        if which_fold(k) > 0; Wk = M_fold{which_fold(k)}.W; else; Wk = M_all.W; end %#ok<PFBNS>
        r = seq_likelihood(w(:, 1), w(:, 2), Wk, ns, k);
        R(k) = struct('L', r.L, 'L_norm', r.L_norm, 'z', r.z, 'pct', r.pct, 'n_spikes', r.n_spikes, ...
            'n_neurons', r.n_neurons);
    end
    sc.name = name;
    sc.cluster_id = ids;
    sc.order_index = (sum(M_all.p, 2) - 0.5) / (n - 1);
    sc.event = ev;                     % index into cand
    sc.which_fold = which_fold;
    sc.L = [R.L]'; sc.L_norm = [R.L_norm]'; sc.z = [R.z]'; sc.pct = [R.pct]';
    sc.n_spikes = [R.n_spikes]'; sc.n_neurons = [R.n_neurons]';
    sc.scored = ~isnan(sc.z);
    sc.forward = sc.pct > 1 - opts.alpha;
    sc.reverse = sc.pct < opts.alpha;
    sc.windows = ev_w;
    k = sc.scored;
    sc.p_binom_forward = 1 - binocdf(sum(sc.forward) - 1, sum(k), opts.alpha);
    sc.p_binom_reverse = 1 - binocdf(sum(sc.reverse) - 1, sum(k), opts.alpha);
    fprintf('%s: %d/%d events scored (>= %d neurons); forward %d (%.1f%%, p = %.2g), reverse %d (%.1f%%, p = %.2g)\n', ...
        name, sum(k), numel(ev), opts.min_neurons, sum(sc.forward), 100 * mean(sc.forward(k)), sc.p_binom_forward, ...
        sum(sc.reverse), 100 * mean(sc.reverse(k)), sc.p_binom_reverse);

    % matched-duration grey control: the cell-ID shuffle is liberal when spikes cluster
    % in time, so every event is also compared with random grey windows of its duration
    rng(opts.seed + p);
    gs = grey_starts(cand.duration(ev), opts.n_control, trials, cand, opts);
    gdur = repelem(cand.duration(ev), opts.n_control);
    ctrl_w = collect_spike_windows(units.spike_times, units.spike_id, ids, gs, gs + gdur);
    zc = nan(numel(ctrl_w), 1);
    Wa = M_all.W;
    nsc = opts.n_shuffle_control;
    parfor j = 1:numel(ctrl_w)
        w = ctrl_w{j};
        if numel(unique(w(:, 2))) < min_n; continue; end
        r = seq_likelihood(w(:, 1), w(:, 2), Wa, nsc, j);
        zc(j) = r.z;
    end
    sc.control_z = zc;                                     % n_control per event, event-major
    sc.control_z95 = prctile(zc, 100 * (1 - opts.alpha));
    sc.control_z05 = prctile(zc, 100 * opts.alpha);
    sc.control_frac_forward = mean(zc(~isnan(zc)) > norminv(1 - opts.alpha));
    sc.forward_cal = sc.forward & sc.z > sc.control_z95;   % also beyond the grey windows' tail
    sc.reverse_cal = sc.reverse & sc.z < sc.control_z05;
    sc.in_stimulus = phase_label(cand.centre(ev), trials, opts.grey_margin);  % A onset to D offset + margin
    g = k & ~sc.in_stimulus;
    sc.p_binom_forward_cal_grey = 1 - binocdf(sum(sc.forward_cal & g) - 1, sum(g), opts.alpha);
    sc.p_binom_reverse_cal_grey = 1 - binocdf(sum(sc.reverse_cal & g) - 1, sum(g), opts.alpha);
    sc.p_binom_forward_cal = 1 - binocdf(sum(sc.forward_cal) - 1, sum(k), opts.alpha);
    sc.p_binom_reverse_cal = 1 - binocdf(sum(sc.reverse_cal) - 1, sum(k), opts.alpha);
    fprintf(['   grey control (%d windows scored): z %g-%g%% = [%.2f %.2f]; calibrated forward %d (p = %.2g), ' ...
        'reverse %d (p = %.2g)\n'], sum(~isnan(zc)), 100 * opts.alpha, 100 * (1 - opts.alpha), sc.control_z05, ...
        sc.control_z95, sum(sc.forward_cal), ...
        sc.p_binom_forward_cal, sum(sc.reverse_cal), sc.p_binom_reverse_cal);
    fprintf('   events in grey only (n = %d): calibrated forward %d (p = %.2g), reverse %d (p = %.2g)\n', sum(g), ...
        sum(sc.forward_cal & g), sc.p_binom_forward_cal_grey, sum(sc.reverse_cal & g), sc.p_binom_reverse_cal_grey);
    scores(p) = sc; %#ok<AGROW>
end

%% 4. timing relative to A onset
rel = @(t) rel_to_onset(t, trials.seq_onset, opts.peri);
timing.ripples = rel(rip.peaktimes(:));
timing.candidates = rel(cand.centre(cand.is_ripple));
timing.candidates_all = rel(cand.centre);
for p = 1:numel(scores)
    tc = cand.centre(scores(p).event);
    timing.(scores(p).name).scored = rel(tc(scores(p).scored));
    timing.(scores(p).name).forward = rel(tc(scores(p).forward));
    timing.(scores(p).name).reverse = rel(tc(scores(p).reverse));
end

%% 5. ripples inside candidate events (one MUA event can hold several ripples)
cand.n_ripples = arrayfun(@(k) sum(rip.peaktimes >= cand.onset(k) & rip.peaktimes <= cand.offset(k)), ...
    (1:numel(cand.onset))');
rip.in_candidate = zeros(numel(rip.peaktimes), 1);       % index of the ripple candidate holding it (0 = none)
rip.in_replay = struct();                                 % per population: inside a significant event
for k = find(cand.is_ripple)'
    rip.in_candidate(rip.peaktimes >= cand.onset(k) & rip.peaktimes <= cand.offset(k)) = k;
end
for p = 1:numel(scores)
    sig_ev = scores(p).event(sig_mask(scores(p), opts));
    rip.in_replay.(scores(p).name) = ismember(rip.in_candidate, sig_ev) & rip.in_candidate > 0;
end
fprintf('Ripples inside ripple-candidate events: %d of %d; ripples per candidate: %s (events with 0/1/2/3+)\n', ...
    sum(rip.in_candidate > 0), numel(rip.peaktimes), mat2str(histcounts(cand.n_ripples(cand.is_ripple), [0 1 2 3 inf])));

%% 6. peri-event population PSTHs (CA1 and V1) around ripple candidates
psth = peri_event_psth(units, scores, cand, trials, opts);

%% figures
figs = [plot_summary(rip, cand, scores, timing, trials, opts, tag), ...
        plot_examples(scores(1), cand, ca1, tvec, LFP.fs, trials, opts, tag)];
names = {'replay_summary', 'replay_examples_CA1'};
if numel(scores) > 1
    figs(end + 1) = plot_examples(scores(2), cand, ca1, tvec, LFP.fs, trials, opts, tag);
    names{end + 1} = sprintf('replay_examples_%s', scores(2).name);
end
figs(end + 1) = plot_peri_psth(psth, scores, opts, tag);
names{end + 1} = 'replay_peri_event_psth';
base = sprintf('replay_%s', stimulus_name);
names = strcat(base, '_', strrep(names, 'replay_', ''));

results.subject = subject;
results.date = date;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.ripples = rip;
results.candidates = cand;
for p = 1:numel(scores); scores(p).windows = []; end
results.scores = scores;
results.timing = timing;
results.psth = psth;
results.ripple_channel = layout.ripple_channel;
if ~opts.save; return; end
safe_save(fullfile(out_dir, [base '_results.mat']), struct('results', results), cfg);
save_vector_figures(figs, names, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function opts = fill_defaults(opts)
def.populations = {'CA1_all', 'V1_good'};   % first = CA1 (MUA population)
def.ripple_peak_min = 5;
def.n_shuffle = 1000;
def.alpha = 0.025;                  % per tail (forward / reverse): 5% two-tailed. The first run used 0.05.
def.psth_window = [-0.5 0.5];       % s around the event onset
def.psth_bin = 0.01;
def.psth_sigma = 0.01;              % s, Gaussian smoothing for display
def.psth_significance = 'shuffle';  % 'shuffle' or 'calibrated' (shuffle + beyond grey 5/95%)
def.min_neurons = 5;
def.peri = [-0.5 1];
def.peri_bin = 0.05;
def.n_examples = 12;
def.n_control = 20;                 % matched-duration grey windows per event
def.n_shuffle_control = 200;
def.grey_margin = 0.3;              % s after D offset before grey control windows
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
end


function m = sig_mask(sc, opts)
% significant replay in either direction (per-tail alpha), optionally grey-calibrated
if strcmp(opts.psth_significance, 'calibrated')
    m = sc.forward_cal | sc.reverse_cal;
else
    m = sc.forward | sc.reverse;
end
end


function psth = peri_event_psth(units, scores, cand, trials, opts)
% Mean firing rate per neuron (Hz) of each population around the onset of every ripple
% candidate, and around matched random grey times (same number, no candidate events).
ev = scores(1).event;
t0 = cand.onset(ev);
rng(opts.seed + 100);
tg = grey_starts(0.01 * ones(numel(ev), 1), 1, trials, cand, opts);
edges = opts.psth_window(1):opts.psth_bin:opts.psth_window(2) + 1e-9;
psth.t = edges(1:end - 1) + opts.psth_bin / 2;
psth.event = ev;
psth.in_stimulus = scores(1).in_stimulus;
psth.event_duration = cand.duration(ev);
for p = 1:numel(scores)
    psth.sig.(scores(p).name) = sig_mask(scores(p), opts);
    ids = scores(p).cluster_id;
    psth.pop(p).name = scores(p).name;
    psth.pop(p).rate = rates(t0);
    psth.pop(p).control = rates(tg);
end

    function R = rates(t)
        w = collect_spike_windows(units.spike_times, units.spike_id, ids, t + opts.psth_window(1), ...
            t + opts.psth_window(2));
        R = zeros(numel(t), numel(edges) - 1);
        for k = 1:numel(t)
            if ~isempty(w{k}); R(k, :) = histcounts(w{k}(:, 1) + opts.psth_window(1), edges); end
        end
        R = R / (numel(ids) * opts.psth_bin);
    end
end


function gs = grey_starts(dur, n_per, trials, cand, opts)
% random window starts in the grey (>= grey_margin after D offset, before the next A),
% not overlapping any MUA candidate event; n_per windows per duration in dur
g0 = trials.seq_offset(1:end - 1) + opts.grey_margin;
g1 = trials.seq_onset(2:end);
gs = zeros(numel(dur) * n_per, 1);
i = 0;
for k = 1:numel(dur)
    span = max(g1 - g0 - dur(k), 0);
    for r = 1:n_per
        for attempt = 1:100
            j = randsample(numel(span), 1, true, span / sum(span));
            s = g0(j) + rand * span(j);
            if ~any(cand.onset < s + dur(k) & cand.offset > s); break; end
        end
        i = i + 1;
        gs(i) = s;
    end
end
end


function s = struct_rows(s)
for f = fieldnames(s)'
    if isnumeric(s.(f{1})) && ~isscalar(s.(f{1})); s.(f{1}) = s.(f{1})(:); end
end
end


function r = rel_to_onset(t, onsets, peri)
% time of each event relative to the A onset it falls near (within peri), NaN otherwise
r = nan(size(t));
for k = 1:numel(t)
    d = t(k) - onsets;
    d = d(d >= peri(1) & d < peri(2));
    if ~isempty(d); r(k) = d(1); end
end
r = r(~isnan(r));
end


%% ------------------------------------------------------------------------
%% Plotting
function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function shade_elements(ax, trials, yl)
c = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.290 0.227 0.655];
for e = 1:numel(trials.element_names)
    x0 = 1000 * trials.element_start(e);
    x1 = x0 + 1000 * trials.element_duration;
    patch(ax, [x0 x1 x1 x0], yl([1 1 2 2]), c(e, :), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax, (x0 + x1) / 2, yl(2), trials.element_names{e}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
end


function f = plot_summary(rip, cand, scores, timing, trials, opts, tag)
n_seq = numel(trials.seq_onset);
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.55 0.55 0.55];
f = figure('Color', 'w', 'Position', [40 40 1600 800]);
tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - awake ripples, MUA candidate events (z >= 3, ripple peak z >= %g) and sequence-likelihood ' ...
    'replay of the ABCD order (%d cell-ID shuffles)'], tag, opts.ripple_peak_min, opts.n_shuffle), 'FontSize', 11);
edges = opts.peri(1):opts.peri_bin:opts.peri(2);
ctr = 1000 * (edges(1:end - 1) + opts.peri_bin / 2);
rate = @(x) histcounts(x, edges) / (n_seq * opts.peri_bin);

% row 1: timing (ripples / candidates, then significant events per population)
ax = nexttile(tl, [1 2]); hold(ax, 'on');
r1 = rate(timing.ripples); r2 = rate(timing.candidates);
yl = [0 1.25 * max([r1 r2 0.05])];
shade_elements(ax, trials, yl);
stairs(ax, 1000 * edges, [r1 r1(end)], '-', 'Color', grey, 'LineWidth', 2, 'DisplayName', ...
    sprintf('ripples (n = %d in window)', numel(timing.ripples)));
stairs(ax, 1000 * edges, [r2 r2(end)], '-', 'Color', blue, 'LineWidth', 2, 'DisplayName', ...
    sprintf('candidate events (n = %d)', numel(timing.candidates)));
ylim(ax, yl); xlim(ax, 1000 * opts.peri);
xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'events / s (per sequence)');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, sprintf('ripples (%d total) and ripple candidate events (%d of %d MUA events) around A onset', ...
    numel(rip.peaktimes), sum(cand.is_ripple), numel(cand.onset)), 'FontWeight', 'normal');
style(ax);

for p = 1:numel(scores)
    ax = nexttile(tl); hold(ax, 'on');
    T = timing.(scores(p).name);
    rf = rate(T.forward); rr = rate(T.reverse); rs = rate(T.scored);
    yl = [0 1.25 * max([rs 0.05])];
    shade_elements(ax, trials, yl);
    stairs(ax, 1000 * edges, [rs rs(end)], '-', 'Color', grey, 'LineWidth', 1.5, 'DisplayName', 'scored candidates');
    stairs(ax, 1000 * edges, [rf rf(end)], '-', 'Color', blue, 'LineWidth', 2, 'DisplayName', ...
        sprintf('forward (n = %d)', numel(T.forward)));
    stairs(ax, 1000 * edges, [rr rr(end)], '-', 'Color', orange, 'LineWidth', 2, 'DisplayName', ...
        sprintf('reverse (n = %d)', numel(T.reverse)));
    ylim(ax, yl); xlim(ax, 1000 * opts.peri);
    xlabel(ax, 'time from A onset (ms)'); ylabel(ax, 'events / s (per sequence)');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s order: significant events (shuffle)', strrep(scores(p).name, '_', ' ')), 'FontWeight', 'normal');
    style(ax);
end

% row 2: per population z distribution (+ grey control) and fractions by phase
for p = 1:numel(scores)
    sc = scores(p);
    k = sc.scored;
    ax = nexttile(tl); hold(ax, 'on');
    zc = sc.control_z(~isnan(sc.control_z));
    e = linspace(min([sc.z(k); zc; -3]), max([sc.z(k); zc; 3]), 35);
    histogram(ax, zc, e, 'Normalization', 'probability', 'FaceColor', grey, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
        'DisplayName', sprintf('grey windows, matched duration (n = %d)', numel(zc)));
    histogram(ax, sc.z(k), e, 'Normalization', 'probability', 'FaceColor', blue, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.6, 'DisplayName', sprintf('candidate events (n = %d)', sum(k)));
    xline(ax, norminv([opts.alpha 1 - opts.alpha]), '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, [sc.control_z05 sc.control_z95], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel(ax, 'L_{seq} z vs cell-ID shuffle'); ylabel(ax, 'fraction');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s (dashed: shuffle 5/95%%, solid: grey 5/95%%)', strrep(sc.name, '_', ' ')), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    tc = cand.centre(sc.event);
    ph = phase_label(tc, trials, opts.grey_margin);
    groups = {'all', sprintf('stimulus (A to D+%.1fs)', opts.grey_margin), 'grey'};
    sel = {k, k & ph == 1, k & ph == 0};
    vals = zeros(numel(groups), 4); nn = zeros(numel(groups), 1);
    for g = 1:numel(groups)
        nn(g) = sum(sel{g});
        vals(g, :) = 100 * [mean(sc.forward(sel{g})) mean(sc.reverse(sel{g})) ...
            mean(sc.forward_cal(sel{g})) mean(sc.reverse_cal(sel{g}))];
    end
    b = bar(ax, vals, 'BarWidth', 0.85, 'EdgeColor', 'none');
    b(1).FaceColor = blue; b(2).FaceColor = orange;
    b(3).FaceColor = blue * 0.55; b(4).FaceColor = orange * 0.6;
    yline(ax, 100 * opts.alpha, ':', 'chance', 'FontSize', 7, 'HandleVisibility', 'off');
    set(ax, 'XTick', 1:numel(groups), 'XTickLabel', strcat(groups, " (n=", string(nn'), ")"));
    ylabel(ax, '% of scored events significant');
    legend(ax, {'forward (shuffle)', 'reverse (shuffle)', 'forward (+grey cal.)', 'reverse (+grey cal.)'}, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 6.5);
    title(ax, sprintf('calibrated, grey events only: forward p = %.2g, reverse p = %.2g', ...
        sc.p_binom_forward_cal_grey, sc.p_binom_reverse_cal_grey), 'FontWeight', 'normal');
    style(ax);
end
end


function f = plot_peri_psth(psth, scores, opts, tag)
% rows: populations; columns: all ripple candidates / candidates in the grey screen only
names = {scores.name};
groups = {'all candidates', [strrep(names{1}, '_', ' ') ' significant'], ...
    [strrep(names{min(2, end)}, '_', ' ') ' significant']};
gcol = [0.15 0.15 0.15; 0.922 0.408 0.204; 0.165 0.471 0.839];
f = figure('Color', 'w', 'Position', [40 40 1300 330 * numel(psth.pop)]);
tl = tiledlayout(f, numel(psth.pop), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - population rate around ripple-candidate onsets (mean per neuron \\pm SEM over events; ' ...
    '%s significance, \\alpha = %g per tail); dashed = matched grey times'], tag, opts.psth_significance, opts.alpha), ...
    'FontSize', 10.5);
k = max(1, round(opts.psth_sigma / opts.psth_bin));
g = exp(-0.5 * ((-3 * k:3 * k) / k).^2); g = g / sum(g);
sm = @(x) conv(x, g, 'same');
t = 1000 * psth.t;
for p = 1:numel(psth.pop)
    for c = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        base = true(numel(psth.event), 1);
        if c == 2; base = ~psth.in_stimulus; end
        masks = {base, base & psth.sig.(names{1}), base & psth.sig.(names{min(2, end)})};
        R = psth.pop(p).rate;
        Cg = psth.pop(p).control(base, :);
        if ~isempty(Cg)
            plot(ax, t, sm(mean(Cg, 1)), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'DisplayName', ...
                sprintf('grey control (n = %d)', size(Cg, 1)));
        end
        for gi = 1:3
            m = masks{gi};
            if sum(m) == 0; continue; end
            mu = sm(mean(R(m, :), 1));
            se = sm(std(R(m, :), 0, 1) / sqrt(sum(m)));
            fill(ax, [t fliplr(t)], [mu - se fliplr(mu + se)], gcol(gi, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
            plot(ax, t, mu, '-', 'Color', gcol(gi, :), 'LineWidth', 1.8, 'DisplayName', ...
                sprintf('%s (n = %d)', groups{gi}, sum(m)));
        end
        xline(ax, 0, 'k:', 'HandleVisibility', 'off');
        xline(ax, 1000 * median(psth.event_duration), 'k:', 'median offset', 'FontSize', 7, 'HandleVisibility', 'off');
        xlim(ax, t([1 end]));
        xlabel(ax, 'time from candidate onset (ms)'); ylabel(ax, 'mean rate per neuron (Hz)');
        ttl = {'all ripple candidates', 'grey-screen candidates only'};
        title(ax, sprintf('%s population: %s', strrep(psth.pop(p).name, '_', ' '), ttl{c}), 'FontWeight', 'normal');
        legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
        set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
    end
end
end


function ph = phase_label(t, trials, margin)
% 1 = stimulus-driven period: A onset to D offset + margin (offset responses included);
% 0 = grey
ph = zeros(size(t));
for k = 1:numel(t)
    ph(k) = any(t(k) >= trials.seq_onset & t(k) < trials.seq_offset + margin);
end
end


function f = plot_examples(sc, cand, lfp, tvec, fs, trials, opts, tag)
% best forward and reverse events: ripple-band LFP above, spikes sorted by the
% template order below (blue = early in the ABCD response, orange = late)
[~, of] = sort(sc.z .* sc.forward, 'descend');
[~, orv] = sort(sc.z .* sc.reverse, 'ascend');
nf = min(sum(sc.forward), ceil(2 * opts.n_examples / 3));
nr = min(sum(sc.reverse), opts.n_examples - nf);
pick = [of(1:nf); orv(1:nr)];
n_ex = numel(pick);
n_col = 4;
n_row = max(1, ceil(n_ex / n_col));
f = figure('Color', 'w', 'Position', [40 40 1600 330 * n_row]);
tl = tiledlayout(f, 3 * n_row, n_col, 'TileSpacing', 'tight', 'Padding', 'compact');
title(tl, sprintf(['%s - %s: strongest forward / reverse replay events (top: 125-300 Hz LFP; bottom: spikes, ' ...
    'neurons sorted by order in the ABCD response, early at bottom)'], tag, strrep(sc.name, '_', ' ')), 'FontSize', 11);
[b, a] = butter(3, [125 300] / (fs / 2), 'bandpass');
[~, rank_order] = sort(sc.order_index);
rank = zeros(numel(rank_order), 1); rank(rank_order) = 1:numel(rank_order);
cmap = interp1([0 0.5 1], [0.165 0.471 0.839; 0.6 0.6 0.6; 0.922 0.408 0.204], linspace(0, 1, 256));
if n_ex == 0
    ax = nexttile(tl); text(ax, 0.5, 0.5, 'no significant events', 'HorizontalAlignment', 'center'); axis(ax, 'off');
    return
end
for i = 1:n_ex
    k = pick(i);
    e = sc.event(k);
    t0 = cand.onset(e); t1 = cand.offset(e);
    pad = 0.05;
    r = floor((i - 1) / n_col); c = mod(i - 1, n_col) + 1;
    ax = nexttile(tl, (3 * r) * n_col + c);
    idx = tvec >= t0 - pad & tvec <= t1 + pad;
    plot(ax, 1000 * (tvec(idx) - t0), filtfilt(b, a, lfp(idx)), 'k', 'LineWidth', 0.5);
    axis(ax, 'tight'); set(ax, 'XTick', []); ylabel(ax, '\muV');
    d = t0 - trials.seq_onset; d = d(d > -1.6 & d < 2);
    when = 'grey';
    if ~isempty(d); [~, j] = min(abs(d)); when = sprintf('%+.0f ms re A', 1000 * d(j)); end
    dirn = 'forward'; if sc.reverse(k); dirn = 'reverse'; end
    if sc.forward_cal(k) || sc.reverse_cal(k); dirn = [dirn ' (grey-cal.)']; end
    title(ax, sprintf('%s z = %.1f (pct %.3f), %s, ripple z %.1f', dirn, sc.z(k), sc.pct(k), when, ...
        cand.ripple_peak(e)), 'FontWeight', 'normal', 'FontSize', 7);
    style(ax);
    ax = nexttile(tl, (3 * r + 1) * n_col + c, [2 1]); hold(ax, 'on');
    w = sc.windows{k};
    if ~isempty(w)
        scatter(ax, 1000 * w(:, 1), rank(w(:, 2)), 14, sc.order_index(w(:, 2)), 'filled');
    end
    colormap(ax, cmap); clim(ax, [min(sc.order_index) max(sc.order_index)]);
    xlim(ax, 1000 * [-pad t1 - t0 + pad]); ylim(ax, [0 numel(rank) + 1]);
    xlabel(ax, 'ms from event onset'); ylabel(ax, 'neuron (template order)');
    style(ax);
end
end
