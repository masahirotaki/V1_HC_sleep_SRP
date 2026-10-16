function results = run_sleep_replay(subject, date, sleep_name, stimulus_name, cfg, opts)
% RUN_SLEEP_REPLAY  Sequence-likelihood (Lseq, Huh et al. 2026) replay of the awake ABCD
% response order in a sleep recording, for ripples, ripple clusters and MUA bursts.
%
%   results = run_sleep_replay('M00069', '20250929', 'Sleep_Box_2', 'GAVNIK_ABCD_1', project_config())
%
% Needs run_abcd_sequence_likelihood on STIMULUS_NAME (units, window) and
% run_sleep_detection on SLEEP_NAME (behavioural_state.mat, detection/extracted_ripple_events.mat).
% Cluster ids come from one sort of the whole day, so ids match across recordings.
%
% Templates: ordered firing probability matrices (seq_order_matrix) from ALL stimulus
%   sequences, per population (opts.populations: CA1_all, V1_good).
% Candidate events (scored when their onset is in opts.states):
%   ripple   every ripple, window [onset, max(offset, onset + 100 ms)]; type solo /
%            cluster / intermediate after Voroslakos et al. 2026 (classify_ripple_clusters:
%            onset-onset < 180 ms = cluster, >= 500 ms from both neighbours = solo)
%   cluster  each ripple cluster as one window: first onset -> end of its last ripple window
%   burst    CA1 MUA bursts with no ripple (detect_mua_bursts: 1 ms, sigma 5 ms, z over SWS,
%            peak z > 3, edges z = 1, <= 300 ms)
% Significance: Lseq vs the cell-ID shuffle (neuron identities permuted across the event's
%   spikes, opts.n_shuffle), forward pct > 1 - alpha, reverse pct < alpha (alpha per tail).
%   Calibrated: also beyond the 1-alpha / alpha percentile of the z of control windows
%   (opts.n_control per event: same duration, same state, no candidate events), because
%   the cell-ID shuffle is liberal when spikes cluster in time.
% Continuation within clusters (spikes inside the member ripple windows only):
%   cross   Lseq of the pairs whose earlier spike is in an earlier member ripple
%           (> 0: later ripples continue the order; < 0: they restart it), z vs cell-ID shuffle
%   shift   mean template order index of the spikes of each member vs member position
%           (slope; > 0: later ripples carry later parts of the sequence), z vs the same shuffle
%   Raw Lseq and its z grow with the number of spikes, so merged vs single z is not compared.
%
% Outputs: <analysis>/<subject>/<date>/<sleep_name>/replay/replay_<sleep>_template_<stimulus>_*

if nargin < 5 || isempty(cfg); cfg = project_config(); end
if nargin < 6; opts = struct(); end
def.populations = {'CA1_all', 'V1_good'};
def.n_shuffle = 1000;
def.alpha = 0.025;
def.min_neurons = 5;
def.n_control = 10;
def.n_shuffle_control = 200;
def.states = {'SWS', 'quietWake'};
def.cluster_iei = 0.18;
def.solo_iei = 0.5;
def.min_window = 0.1;
def.mua_sigma = 0.005;
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
tag = sprintf('%s %s %s (template: %s)', subject, date, strrep(sleep_name, '_', ' '), strrep(stimulus_name, '_', ' '));

%% templates from the stimulus recording
spaths  = get_session_paths(subject, date, stimulus_name, cfg);
sunits  = load_session_clusters(spaths, get_probe_layout(spaths, cfg), cfg);
strials = get_abcd_trials(spaths);
S = load(fullfile(spaths.analysis_dir, 'Decoding', sprintf('decoding_%s_seqlik_results.mat', stimulus_name)));
win = S.results.opts.window;

%% sleep recording: units, states, ripples, MUA bursts
paths  = get_session_paths(subject, date, sleep_name, cfg);
units  = load_session_clusters(paths, get_probe_layout(paths, cfg), cfg);
out_dir = fullfile(paths.analysis_dir, 'replay');
check_output_path(out_dir, cfg);
B = load(fullfile(paths.output_dir, 'behavioural_state.mat'));
states = B.behavioural_state;
R = load(fullfile(paths.detection_dir, 'extracted_ripple_events.mat'));
rip = R.ripples;
C = classify_ripple_clusters(rip.onset(:), rip.offset(:), struct('cluster_iei', opts.cluster_iei, ...
    'solo_iei', opts.solo_iei, 'min_window', opts.min_window));
ca1 = units.cluster_id(units.region == "CA1");
sp = units.spike_times(ismember(units.spike_id, ca1));
t_end = max([sp; C.offset]) + 1;
[Bu, ~] = detect_mua_bursts(sp, t_end, struct('sigma', opts.mua_sigma, 'baseline', states.SWS));
has_rip = arrayfun(@(i) any(C.window(:, 1) < Bu.offset(i) & C.window(:, 2) > Bu.onset(i)), (1:numel(Bu.onset))');

% event table
nr = numel(C.onset); nc = numel(C.clusters.onset); nb = sum(~has_rip);
E = table();
E.kind = [repmat("ripple", nr, 1); repmat("cluster", nc, 1); repmat("burst", nb, 1)];
E.onset = [C.window(:, 1); C.clusters.onset; Bu.onset(~has_rip)];
E.offset = [C.window(:, 2); C.clusters.offset; Bu.offset(~has_rip)];
E.type = [C.type; repmat("cluster", nc, 1); repmat("burst", nb, 1)];
E.ripple_index = [(1:nr)'; zeros(nc + nb, 1)];
E.cluster_id = [C.cluster; (1:nc)'; zeros(nb, 1)];
E.position = [C.position; zeros(nc + nb, 1)];
E.n_ripples = [ones(nr, 1); C.clusters.size; zeros(nb, 1)];
E.state = state_at(E.onset, states);
E.duration = E.offset - E.onset;
use = ismember(E.state, string(opts.states));
fprintf('%s: %d ripples (%s), %d clusters, %d MUA bursts without ripple; scoring %d events in %s\n', sleep_name, nr, ...
    strjoin(arrayfun(@(t) sprintf('%s %d', t, sum(C.type == t)), ["solo" "cluster" "intermediate"], 'UniformOutput', false), ', '), ...
    nc, nb, sum(use), strjoin(opts.states, '/'));

%% scoring
ev = find(use);
cand_iv = [C.window; [Bu.onset Bu.offset]];
for p = 1:numel(opts.populations)
    name = opts.populations{p};
    ids = S.results.pop(strcmp({S.results.pop.name}, name)).cluster_id;
    n = numel(ids);
    tw = collect_spike_windows(sunits.spike_times, sunits.spike_id, ids, strials.seq_onset + win(1), strials.seq_onset + win(2));
    M = seq_order_matrix(tw, n);
    W = M.W;
    oi = (sum(M.p, 2) - 0.5) / (n - 1);
    ew = collect_spike_windows(units.spike_times, units.spike_id, ids, E.onset(ev), E.offset(ev));
    [z, pct, L, Ln, ns, nn] = deal(nan(numel(ev), 1));
    nsh = opts.n_shuffle; min_n = opts.min_neurons;
    parfor k = 1:numel(ev)
        w = ew{k};
        nn(k) = numel(unique(w(:, 2))); ns(k) = size(w, 1);
        if nn(k) < min_n; continue; end
        r = seq_likelihood(w(:, 1), w(:, 2), W, nsh, k);
        z(k) = r.z; pct(k) = r.pct; L(k) = r.L; Ln(k) = r.L_norm;
    end
    % matched-duration controls in the same state, outside candidate events
    rng(opts.seed + p);
    [cs, cstate] = control_starts(E.duration(ev), E.state(ev), states, cand_iv, opts.n_control);
    cw = collect_spike_windows(units.spike_times, units.spike_id, ids, cs, cs + repelem(E.duration(ev), opts.n_control));
    zc = nan(numel(cw), 1); nsc = opts.n_shuffle_control;
    parfor j = 1:numel(cw)
        w = cw{j};
        if numel(unique(w(:, 2))) < min_n; continue; end
        r = seq_likelihood(w(:, 1), w(:, 2), W, nsc, j);
        zc(j) = r.z;
    end
    sc = struct('name', name, 'cluster_id', ids, 'order_index', oi, 'event', ev, 'z', z, 'pct', pct, 'L', L, ...
        'L_norm', Ln, 'n_spikes', ns, 'n_neurons', nn, 'scored', ~isnan(z), 'forward', pct > 1 - opts.alpha, ...
        'reverse', pct < opts.alpha, 'control_z', zc, 'control_state', cstate);
    sc.forward_cal = false(numel(ev), 1); sc.reverse_cal = false(numel(ev), 1);
    sc.control_hi = struct(); sc.control_lo = struct(); sc.control_frac_sig = struct();
    for s = string(opts.states)
        es = E.state(ev) == s; czs = zc(cstate == s & ~isnan(zc));
        if numel(czs) < 20; continue; end
        sc.forward_cal(es) = sc.forward(es) & z(es) > prctile(czs, 100 * (1 - opts.alpha));
        sc.reverse_cal(es) = sc.reverse(es) & z(es) < prctile(czs, 100 * opts.alpha);
        sc.control_hi.(s) = prctile(czs, 100 * (1 - opts.alpha));
        sc.control_lo.(s) = prctile(czs, 100 * opts.alpha);
        sc.control_frac_sig.(s) = mean(czs > norminv(1 - opts.alpha) | czs < norminv(opts.alpha));
    end
    % continuation inside clusters
    sc.continuation = continuation_tests(C, E, ev, units, ids, W, oi, opts);
    for kind = ["ripple" "cluster" "burst"]
        m = sc.scored & E.kind(ev) == kind;
        fprintf('  %s %-8s: %4d scored, forward %3d (%.1f%%), reverse %3d (%.1f%%); calibrated %d / %d\n', name, kind, ...
            sum(m), sum(sc.forward(m)), 100 * mean(sc.forward(m)), sum(sc.reverse(m)), 100 * mean(sc.reverse(m)), ...
            sum(sc.forward_cal(m)), sum(sc.reverse_cal(m)));
    end
    scores(p) = sc; %#ok<AGROW>
end

results.subject = subject;
results.date = date;
results.sleep = sleep_name;
results.stimulus = stimulus_name;
results.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
results.opts = opts;
results.events = E;
results.clusters = C;
results.mua_bursts = Bu;
results.scores = scores;
results.summary = summary_table(scores, E, opts);
disp(results.summary);
figs = [plot_summary(scores, E, opts, tag), plot_continuation(scores, opts, tag)];
if ~opts.save; return; end
base = sprintf('replay_%s_template_%s', sleep_name, stimulus_name);
safe_save(fullfile(out_dir, [base '_results.mat']), struct('results', results), cfg);
csv = fullfile(out_dir, [base '_summary.csv']);
check_output_path(csv, cfg); check_overwrite(csv, cfg);
writetable(results.summary, csv);
save_vector_figures(figs, {[base '_summary'], [base '_continuation']}, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function st = state_at(t, states)
st = repmat("other", numel(t), 1);
for s = ["SWS" "REM" "quietWake" "movement"]
    if ~isfield(states, s); continue; end
    iv = states.(s);
    for k = 1:size(iv, 1); st(t >= iv(k, 1) & t < iv(k, 2)) = s; end
end
end

function [cs, cstate] = control_starts(dur, est, states, cand_iv, n_per)
% n_per random starts per event: same duration, inside the same state, not overlapping
% any candidate interval; NaN when the state has no room (-> empty window, NaN z)
cs = nan(numel(dur) * n_per, 1);
cstate = repelem(est(:), n_per);
i = 0;
for k = 1:numel(dur)
    if est(k) == "other" || ~isfield(states, est(k)); i = i + n_per; continue; end
    iv = states.(est(k));
    len = max(iv(:, 2) - iv(:, 1) - dur(k), 0);
    for r = 1:n_per
        i = i + 1;
        if sum(len) == 0; continue; end
        for attempt = 1:100
            j = randsample(numel(len), 1, true, len / sum(len));
            t = iv(j, 1) + rand * len(j);
            if ~any(cand_iv(:, 1) < t + dur(k) & cand_iv(:, 2) > t); break; end
        end
        cs(i) = t;
    end
end
end

function T = continuation_tests(C, E, ev, units, ids, W, oi, opts)
% clusters with >= 2 member ripples: cross-ripple Lseq term and order-index shift, each
% z-scored against the cell-ID shuffle of the cluster's spikes (member labels and times kept)
cl = find(C.clusters.size >= 2);
cl = cl(ismember(cl, E.cluster_id(ev(E.kind(ev) == "cluster"))));
T = table();
T.cluster_id = cl(:);
[T.n_ripples, T.n_spikes, T.cross, T.cross_z, T.shift, T.shift_z] = deal(nan(numel(cl), 1));
n = numel(ids);
for c = 1:numel(cl)
    mem = C.clusters.members{cl(c)};
    T.n_ripples(c) = numel(mem);
    w = collect_spike_windows(units.spike_times, units.spike_id, ids, C.window(mem, 1), C.window(mem, 2));
    lab = cell2mat(arrayfun(@(m) m * ones(size(w{m}, 1), 1), (1:numel(mem))', 'UniformOutput', false));
    id = cell2mat(cellfun(@(x) x(:, 2), w, 'UniformOutput', false));
    T.n_spikes(c) = numel(id);
    if numel(unique(id)) < opts.min_neurons || numel(unique(lab)) < 2; continue; end
    [x0, s0] = stats(id, lab);
    xs = zeros(opts.n_shuffle, 1); ss = xs;
    for k = 1:opts.n_shuffle
        [xs(k), ss(k)] = stats(id(randperm(numel(id))), lab);
    end
    T.cross(c) = x0; T.cross_z(c) = (x0 - mean(xs)) / std(xs);
    T.shift(c) = s0; T.shift_z(c) = (s0 - mean(ss)) / std(ss);
end

    function [cross, shift] = stats(idv, labv)
        m = max(labv);
        cnt = zeros(n, m);
        for q = 1:m; cnt(:, q) = accumarray(idv(labv == q), 1, [n 1]); end
        prev = [zeros(n, 1), cumsum(cnt(:, 1:end - 1), 2)];   % spikes of earlier members
        cross = sum(sum(cnt .* (W * prev)));                    % sum_b W(id_b, :) * prev
        mo = arrayfun(@(q) mean(oi(idv(labv == q))), 1:m);
        ok = ~isnan(mo);
        if sum(ok) >= 2
            pp = polyfit(find(ok), mo(ok), 1); shift = pp(1);
        else
            shift = NaN;
        end
    end
end

function T = summary_table(scores, E, opts)
rows = {};
for p = 1:numel(scores)
    sc = scores(p);
    st = E.state(sc.event); kd = E.kind(sc.event); ty = E.type(sc.event);
    groups = {"ripple: solo", kd == "ripple" & ty == "solo"; "ripple: in cluster", kd == "ripple" & ty == "cluster"; ...
        "ripple: intermediate", kd == "ripple" & ty == "intermediate"; "cluster (merged)", kd == "cluster"; ...
        "MUA burst, no ripple", kd == "burst"};
    for s = string(opts.states)
        for g = 1:size(groups, 1)
            m = sc.scored & st == s & groups{g, 2};
            n = sum(m);
            rows(end + 1, :) = {sc.name, char(s), char(groups{g, 1}), n, mean(sc.forward(m)), mean(sc.reverse(m)), ...
                mean(sc.forward_cal(m)), mean(sc.reverse_cal(m)), 1 - binocdf(sum(sc.forward(m)) - 1, n, opts.alpha), ...
                1 - binocdf(sum(sc.forward_cal(m)) - 1, n, opts.alpha), 1 - binocdf(sum(sc.reverse_cal(m)) - 1, n, opts.alpha)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'population', 'state', 'event', 'n_scored', 'frac_forward', 'frac_reverse', ...
    'frac_forward_cal', 'frac_reverse_cal', 'p_forward', 'p_forward_cal', 'p_reverse_cal'});
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function f = plot_summary(scores, E, opts, tag)
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.55 0.55 0.55];
np = numel(scores);
f = figure('Color', 'w', 'Position', [40 40 1650 360 * np]);
tl = tiledlayout(f, np, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - Lseq replay of the awake ABCD order; cell-ID shuffle (%d), \\alpha = %g per tail; ' ...
    'calibrated = also beyond matched-duration control windows'], tag, opts.n_shuffle, opts.alpha), 'FontSize', 10.5);
for p = 1:np
    sc = scores(p);
    kd = E.kind(sc.event); ty = E.type(sc.event); st = E.state(sc.event);
    ax = nexttile(tl); hold(ax, 'on');
    zc = sc.control_z(~isnan(sc.control_z));
    e = linspace(-5, 6, 45);
    histogram(ax, zc, e, 'Normalization', 'probability', 'FaceColor', grey, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
        'DisplayName', sprintf('control windows (n = %d)', numel(zc)));
    histogram(ax, sc.z(sc.scored & kd == "ripple"), e, 'Normalization', 'probability', 'FaceColor', blue, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', sprintf('ripples (n = %d)', sum(sc.scored & kd == "ripple")));
    histogram(ax, sc.z(sc.scored & kd == "cluster"), e, 'Normalization', 'probability', 'FaceColor', orange, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', sprintf('merged clusters (n = %d)', sum(sc.scored & kd == "cluster")));
    xline(ax, norminv([opts.alpha 1 - opts.alpha]), '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xlabel(ax, 'L_{seq} z vs cell-ID shuffle'); ylabel(ax, 'fraction');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s template', strrep(sc.name, '_', ' ')), 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    gl = {'solo', 'in cluster', 'intermed.', 'merged cluster', 'burst only'};
    gm = {kd == "ripple" & ty == "solo", kd == "ripple" & ty == "cluster", kd == "ripple" & ty == "intermediate", ...
        kd == "cluster", kd == "burst"};
    sws = st == "SWS";
    vals = zeros(numel(gm), 2); nn = zeros(numel(gm), 1);
    for g = 1:numel(gm)
        m = sc.scored & sws & gm{g};
        nn(g) = sum(m);
        vals(g, :) = 100 * [mean(sc.forward_cal(m)) mean(sc.reverse_cal(m))];
    end
    b = bar(ax, vals, 'BarWidth', 0.85, 'EdgeColor', 'none');
    b(1).FaceColor = blue; b(2).FaceColor = orange;
    yline(ax, 100 * opts.alpha, ':', 'chance per tail', 'FontSize', 7);
    set(ax, 'XTick', 1:numel(gl), 'XTickLabel', strcat(gl, " (", string(nn'), ")"));
    ylabel(ax, '% calibrated significant (SWS)');
    legend(ax, {'forward', 'reverse'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, 'SWS events by type', 'FontWeight', 'normal');
    style(ax);

    ax = nexttile(tl); hold(ax, 'on');
    t = E.onset(sc.event) / 60;
    plot(ax, t(sc.scored), sc.z(sc.scored), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 5);
    plot(ax, t(sc.forward_cal), sc.z(sc.forward_cal), 'o', 'Color', blue, 'MarkerFaceColor', blue, 'MarkerSize', 3);
    plot(ax, t(sc.reverse_cal), sc.z(sc.reverse_cal), 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 3);
    yline(ax, 0, 'k:');
    xlabel(ax, 'time in recording (min)'); ylabel(ax, 'L_{seq} z');
    title(ax, 'events across the session (filled = calibrated significant)', 'FontWeight', 'normal');
    style(ax);
end
end

function f = plot_continuation(scores, opts, tag)
np = numel(scores);
f = figure('Color', 'w', 'Position', [40 40 1100 330 * np]);
tl = tiledlayout(f, np, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - inside ripple clusters: does the next ripple continue the sequence? ' ...
    '(z vs cell-ID shuffle, %d)'], tag, opts.n_shuffle), 'FontSize', 10.5);
for p = 1:np
    T = scores(p).continuation;
    for q = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        if q == 1; z = T.cross_z; lab = 'cross-ripple L_{seq} term (z)';
        else; z = T.shift_z; lab = 'template-position shift across member ripples (z)'; end
        z = z(~isnan(z));
        if ~isempty(z)
            histogram(ax, z, linspace(-4, 4, 25), 'FaceColor', [0.165 0.471 0.839], 'EdgeColor', 'none');
            xline(ax, 0, 'k-');
            pv = NaN; if numel(z) >= 5; pv = signrank(z); end
            title(ax, sprintf('%s: n = %d clusters, mean %.2f, signrank p = %.2g', strrep(scores(p).name, '_', ' '), ...
                numel(z), mean(z), pv), 'FontWeight', 'normal');
        else
            title(ax, sprintf('%s: no scorable clusters', strrep(scores(p).name, '_', ' ')), 'FontWeight', 'normal');
        end
        xlabel(ax, lab); ylabel(ax, 'clusters');
        style(ax);
    end
end
end
