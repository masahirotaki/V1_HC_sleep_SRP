function results = analyze_ripple_clusters(subject, date, sleep_name, cfg, opts)
% ANALYZE_RIPPLE_CLUSTERS  How to define replay windows around sleep ripples: singlets vs
% ripple clusters (doublets, triplets), CA1 MUA after the last ripple, spikes available
% per window, and whether ripples in one cluster carry similar content.
%
%   results = analyze_ripple_clusters('M00069', '20250929', 'Sleep_Box_1', project_config())
%
% Needs run_sleep_detection on SLEEP_NAME (detection/extracted_ripple_events.mat,
% behavioural_state.mat).
%
% Clusters: ripples are chained when the next ripple starts less than opts.cluster_gap
% after the previous one ends (offset -> next onset). Candidate windows compared:
%   'onset+100'    [first onset, first onset + 100 ms]          (paper style, singlets)
%   'offset+100'   [first onset, last offset + 100 ms]          (cluster-aware)
%   'mua'          [first onset, CA1 MUA back below z = 0 after the last offset]
% Content similarity between two ripples (template-free; spikes onset -> offset):
%   population vector: cosine similarity of spike counts per neuron
%   order: Spearman rho of the neurons' mean spike time (re each ripple's onset) over the
%          neurons active in both (>= opts.min_shared), z-scored against shuffled neuron labels
% compared for consecutive ripples within a cluster, neighbouring separate ripples (onset
% gap 0.5-5 s) and random ripple pairs (> 60 s apart), all in SWS.
%
% Outputs: <analysis>/<subject>/<date>/<sleep_name>/replay/replay_<sleep_name>_ripple_clusters.*

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.cluster_gap = 0.1;             % s, offset -> next onset
def.post = 0.1;                    % s added after the (last) ripple
def.mua_bin = 0.001;
def.mua_window = 0.02;             % s, gausswin length for MUA (0.1 = detect_candidate_events_masa)
def.min_neurons = 5;
def.min_shared = 5;
def.n_shuffle = 500;
def.n_random_pairs = 1000;
def.populations = {'CA1', 'V1'};   % CA1: all units; V1: good units
def.seed = 1;
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
rng(opts.seed);

paths  = get_session_paths(subject, date, sleep_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
R = load(fullfile(paths.detection_dir, 'extracted_ripple_events.mat'));
rip = R.ripples;
B = load(fullfile(paths.output_dir, 'behavioural_state.mat'));
states = B.behavioural_state;
tag = sprintf('%s %s %s', subject, date, strrep(sleep_name, '_', ' '));

on = rip.onset(:); off = rip.offset(:); st = string(rip.state(:));
[on, o] = sort(on); off = off(o); st = st(o);
n = numel(on);
dur = off - on;
iri = diff(on);                       % onset -> next onset
gap = on(2:end) - off(1:end - 1);     % offset -> next onset
sws = st == "SWS";
fprintf('%s: %d ripples (%d SWS); duration median %.0f ms; onset-onset median %.2f s\n', sleep_name, n, sum(sws), ...
    1000 * median(dur), median(iri));

%% clusters
cl = cumsum([1; gap >= opts.cluster_gap]);      % cluster id per ripple
nc = max(cl);
csize = accumarray(cl, 1);
c_on = accumarray(cl, on, [], @min);
c_off = accumarray(cl, off, [], @max);
c_state = arrayfun(@(k) st(find(cl == k, 1)), (1:nc)');
within = cl(1:end - 1) == cl(2:end);             % consecutive pair inside one cluster
fprintf('clusters (gap < %.0f ms): %d events; singlets %d, doublets %d, triplets %d, 4+ %d; %.0f%% of ripples in clusters\n', ...
    1000 * opts.cluster_gap, nc, sum(csize == 1), sum(csize == 2), sum(csize == 3), sum(csize >= 4), ...
    100 * mean(csize(cl) > 1));
gaps_grid = (0:0.01:0.5)';
n_events = arrayfun(@(g) 1 + sum(gap >= g), gaps_grid);

%% CA1 MUA (z over SWS)
ca1 = units.cluster_id(units.region == "CA1");
sp = units.spike_times(ismember(units.spike_id, ca1));
T = max([off; sp]) + 1;
edges = 0:opts.mua_bin:T;
tm = edges(1:end - 1) + opts.mua_bin / 2;
cnt = histcounts(sp, edges);
g = gausswin(max(3, round(opts.mua_window / opts.mua_bin)))'; g = g / sum(g);
mua = conv(cnt, g, 'same');
in_sws = false(size(tm));
for k = 1:size(states.SWS, 1); in_sws(tm >= states.SWS(k, 1) & tm < states.SWS(k, 2)) = true; end
mz = (mua - mean(mua(in_sws))) / std(mua(in_sws));
% time from the last ripple offset until MUA z < 0 (searched up to 0.5 s)
t_back = nan(nc, 1);
for k = 1:nc
    i0 = find(tm >= c_off(k), 1);
    i1 = find(mz(i0:min(end, i0 + round(0.5 / opts.mua_bin))) < 0, 1);
    if ~isempty(i1); t_back(k) = (i1 - 1) * opts.mua_bin; end
end
% peri-onset MUA, singlets vs doublets (SWS)
pw = -0.2:opts.mua_bin * 5:0.5;
peri = @(t) cell2mat(arrayfun(@(x) interp1(tm, mz, x + pw, 'linear', NaN), t(:), 'UniformOutput', false));
mua_singlet = peri(c_on(csize == 1 & c_state == "SWS"));
mua_doublet = peri(c_on(csize == 2 & c_state == "SWS"));
mua_multi = peri(c_on(csize >= 3 & c_state == "SWS"));

%% windows and spikes available
win.onset100 = [c_on, c_on + 0.1];
win.offset100 = [c_on, c_off + opts.post];
win.mua = [c_on, c_off + fillmissing(t_back, 'constant', opts.post)];
wn = fieldnames(win);
avail = struct();
for p = 1:numel(opts.populations)
    ids = pop_ids(units, opts.populations{p});
    for w = 1:numel(wn)
        ww = collect_spike_windows(units.spike_times, units.spike_id, ids, win.(wn{w})(:, 1), win.(wn{w})(:, 2));
        avail.(opts.populations{p}).(wn{w}).n_neurons = cellfun(@(x) numel(unique(x(:, 2))), ww);
        avail.(opts.populations{p}).(wn{w}).n_spikes = cellfun(@(x) size(x, 1), ww);
    end
end

%% content similarity between ripples (SWS only)
idx = (1:n)';
pairs.within = [idx(1:end - 1), idx(2:end)];
pairs.within = pairs.within(within & sws(1:end - 1) & sws(2:end), :);
nb = ~within & iri >= 0.5 & iri <= 5 & sws(1:end - 1) & sws(2:end);
pairs.neighbour = [idx(1:end - 1), idx(2:end)];
pairs.neighbour = pairs.neighbour(nb, :);
ss = find(sws);
a = ss(randi(numel(ss), opts.n_random_pairs * 3, 1)); b = ss(randi(numel(ss), opts.n_random_pairs * 3, 1));
far = abs(on(a) - on(b)) > 60;
pairs.random = [a(far), b(far)];
pairs.random = pairs.random(1:min(end, opts.n_random_pairs), :);
pn = fieldnames(pairs);
sim = struct();
for p = 1:numel(opts.populations)
    ids = pop_ids(units, opts.populations{p});
    rw = collect_spike_windows(units.spike_times, units.spike_id, ids, on, off);
    nid = numel(ids);
    V = zeros(n, nid); MT = nan(n, nid);
    for k = 1:n
        if isempty(rw{k}); continue; end
        V(k, :) = accumarray(rw{k}(:, 2), 1, [nid 1])';
        mt = accumarray(rw{k}(:, 2), rw{k}(:, 1), [nid 1], @mean, NaN)';
        MT(k, :) = mt;
    end
    for q = 1:numel(pn)
        P = pairs.(pn{q});
        cs = nan(size(P, 1), 1); rz = cs; rr = cs;
        for k = 1:size(P, 1)
            v1 = V(P(k, 1), :); v2 = V(P(k, 2), :);
            if any(v1) && any(v2); cs(k) = (v1 * v2') / (norm(v1) * norm(v2)); end
            sh = ~isnan(MT(P(k, 1), :)) & ~isnan(MT(P(k, 2), :));
            if sum(sh) >= opts.min_shared
                x = MT(P(k, 1), sh); y = MT(P(k, 2), sh);
                rr(k) = corr(x', y', 'type', 'Spearman');
                null = zeros(opts.n_shuffle, 1);
                for s = 1:opts.n_shuffle; null(s) = corr(x', y(randperm(numel(y)))', 'type', 'Spearman'); end
                rz(k) = (rr(k) - mean(null)) / std(null);
            end
        end
        sim.(opts.populations{p}).(pn{q}) = struct('cosine', cs, 'rho', rr, 'rho_z', rz);
    end
    S = sim.(opts.populations{p});
    fprintf('%s similarity (SWS): cosine within %.2f / neighbour %.2f / random %.2f; order z within %.2f (n=%d) / neighbour %.2f (n=%d) / random %.2f (n=%d)\n', ...
        opts.populations{p}, median(S.within.cosine, 'omitnan'), median(S.neighbour.cosine, 'omitnan'), ...
        median(S.random.cosine, 'omitnan'), mean(S.within.rho_z, 'omitnan'), sum(~isnan(S.within.rho_z)), ...
        mean(S.neighbour.rho_z, 'omitnan'), sum(~isnan(S.neighbour.rho_z)), mean(S.random.rho_z, 'omitnan'), ...
        sum(~isnan(S.random.rho_z)));
end

results = struct('subject', subject, 'date', date, 'sleep', sleep_name, 'opts', opts, 'onset', on, 'offset', off, ...
    'state', st, 'cluster', cl, 'cluster_size', csize, 'cluster_onset', c_on, 'cluster_offset', c_off, ...
    'cluster_state', c_state, 'mua_back_below0', t_back, 'windows', win, 'available', avail, 'pairs', pairs, ...
    'similarity', sim, 'gap_grid', gaps_grid, 'n_events_vs_gap', n_events);
fig = plot_clusters(results, iri, gap, dur, sws, pw, mua_singlet, mua_doublet, mua_multi, tag, opts);
if ~opts.save; return; end
out_dir = fullfile(paths.analysis_dir, 'replay');
base = sprintf('replay_%s_ripple_clusters', sleep_name);
safe_save(fullfile(out_dir, [base '.mat']), struct('results', results), cfg);
save_vector_figures(fig, {base}, out_dir, cfg);
end


%% ------------------------------------------------------------------------
function ids = pop_ids(units, name)
if strcmp(name, 'CA1')
    ids = units.cluster_id(units.region == "CA1");
else
    ids = units.cluster_id(units.region == "V1" & units.is_good);
end
end

function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function f = plot_clusters(r, iri, gap, dur, sws, pw, ms, md, mm, tag, opts)
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; green = [0.106 0.686 0.478]; grey = [0.55 0.55 0.55];
f = figure('Color', 'w', 'Position', [30 30 1750 1000]);
tl = tiledlayout(f, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s - ripple singlets vs clusters (chained if offset \\rightarrow next onset < %.0f ms) and candidate replay windows', ...
    tag, 1000 * opts.cluster_gap), 'FontSize', 11);
s2 = sws(1:end - 1) & sws(2:end);

ax = nexttile(tl); hold(ax, 'on');
histogram(ax, 1000 * dur(sws), 0:10:300, 'FaceColor', blue, 'EdgeColor', 'none');
xline(ax, 100, 'k--', '100 ms', 'FontSize', 7);
xlabel(ax, 'ripple duration (ms)'); ylabel(ax, 'ripples (SWS)');
title(ax, sprintf('duration: median %.0f ms, %.0f%% > 100 ms', 1000 * median(dur(sws)), 100 * mean(dur(sws) > 0.1)), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
e = logspace(log10(0.02), log10(100), 60);
histogram(ax, iri(s2), e, 'FaceColor', blue, 'EdgeColor', 'none');
set(ax, 'XScale', 'log');
xline(ax, [0.1 0.2 1], 'k:');
xlabel(ax, 'onset \rightarrow next onset (s, log)'); ylabel(ax, 'ripple pairs (SWS)');
title(ax, sprintf('inter-ripple interval: %.0f%% < 200 ms, %.0f%% < 1 s', 100 * mean(iri(s2) < 0.2), 100 * mean(iri(s2) < 1)), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
e2 = logspace(log10(0.001), log10(100), 60);
histogram(ax, max(gap(s2), 0.001), e2, 'FaceColor', green, 'EdgeColor', 'none');
set(ax, 'XScale', 'log');
xline(ax, opts.cluster_gap, 'k--', sprintf('%.0f ms', 1000 * opts.cluster_gap), 'FontSize', 7);
xlabel(ax, 'offset \rightarrow next onset gap (s, log)'); ylabel(ax, 'ripple pairs (SWS)');
title(ax, 'gap between ripples (sets the chaining)', 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, 1000 * r.gap_grid, r.n_events_vs_gap, '-', 'Color', blue, 'LineWidth', 2);
xline(ax, 1000 * opts.cluster_gap, 'k--');
xlabel(ax, 'chaining gap threshold (ms)'); ylabel(ax, 'number of events');
title(ax, sprintf('%d ripples \\rightarrow %d events at %.0f ms', numel(r.onset), max(r.cluster), 1000 * opts.cluster_gap), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

% row 2
ax = nexttile(tl); hold(ax, 'on');
cs = r.cluster_size(r.cluster_state == "SWS");
cnts = [sum(cs == 1) sum(cs == 2) sum(cs == 3) sum(cs >= 4)];
b = bar(ax, cnts, 'FaceColor', blue, 'EdgeColor', 'none');
text(ax, b.XEndPoints, cnts, string(cnts), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
set(ax, 'XTick', 1:4, 'XTickLabel', {'singlet', 'doublet', 'triplet', '4+'});
ylabel(ax, 'events (SWS)');
title(ax, sprintf('event types: %.0f%% of SWS ripples are in clusters', ...
    100 * mean(r.cluster_size(r.cluster(r.state == "SWS")) > 1)), 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
w = r.cluster(1:end - 1) == r.cluster(2:end) & s2;
histogram(ax, 1000 * iri(w), 0:10:400, 'FaceColor', orange, 'EdgeColor', 'none');
xlabel(ax, 'onset \rightarrow next onset inside a cluster (ms)'); ylabel(ax, 'ripple pairs');
title(ax, sprintf('within-cluster interval: median %.0f ms (n = %d)', 1000 * median(iri(w)), sum(w)), ...
    'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
cols = {blue, orange, green}; lab = {'singlet', 'doublet', 'triplet+'}; M = {ms, md, mm};
for k = 1:3
    if isempty(M{k}); continue; end
    mu = mean(M{k}, 1, 'omitnan'); se = std(M{k}, 0, 1, 'omitnan') / sqrt(size(M{k}, 1));
    fill(ax, 1000 * [pw fliplr(pw)], [mu - se fliplr(mu + se)], cols{k}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    plot(ax, 1000 * pw, mu, '-', 'Color', cols{k}, 'LineWidth', 2, 'DisplayName', sprintf('%s (n = %d)', lab{k}, size(M{k}, 1)));
end
xline(ax, [0 100], 'k:', 'HandleVisibility', 'off'); yline(ax, 0, 'k-', 'HandleVisibility', 'off');
xlabel(ax, 'time from first ripple onset (ms)'); ylabel(ax, sprintf('CA1 MUA z (%.0f ms window, SWS baseline)', 1000 * opts.mua_window));
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, 'CA1 MUA around events', 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

ax = nexttile(tl); hold(ax, 'on');
tb = r.mua_back_below0(r.cluster_state == "SWS");
histogram(ax, 1000 * tb, 0:10:500, 'FaceColor', grey, 'EdgeColor', 'none');
xline(ax, 1000 * opts.post, 'k--', sprintf('+%.0f ms', 1000 * opts.post), 'FontSize', 7);
xlabel(ax, 'last ripple offset \rightarrow MUA z < 0 (ms)'); ylabel(ax, 'events (SWS)');
title(ax, sprintf('MUA returns to baseline: median %.0f ms after offset (%d never within 0.5 s)', 1000 * median(tb, 'omitnan'), ...
    sum(isnan(tb))), 'FontWeight', 'normal', 'FontSize', 8.5);
style(ax);

% row 3: availability per window, and similarity
wn = fieldnames(r.windows); wl = {'onset + 100 ms', 'offset + 100 ms', 'until MUA < 0'};
ws = r.cluster_state == "SWS";
for p = 1:numel(opts.populations)
    ax = nexttile(tl); hold(ax, 'on');
    A = r.available.(opts.populations{p});
    for w = 1:numel(wn)
        nn = A.(wn{w}).n_neurons(ws);
        [cdf, x] = ecdf(nn);
        stairs(ax, x, 1 - cdf, 'LineWidth', 2, 'DisplayName', sprintf('%s: %.0f%% \\geq %d, med dur %.0f ms', wl{w}, ...
            100 * mean(nn >= opts.min_neurons), opts.min_neurons, 1000 * median(diff(r.windows.(wn{w})(ws, :), [], 2))));
    end
    xline(ax, opts.min_neurons, 'k:', 'HandleVisibility', 'off');
    xlabel(ax, sprintf('active %s neurons in the window', opts.populations{p})); ylabel(ax, 'fraction of SWS events \geq x');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, sprintf('%s neurons available per window', opts.populations{p}), 'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
end
for p = 1:numel(opts.populations)
    ax = nexttile(tl); hold(ax, 'on');
    S = r.similarity.(opts.populations{p});
    groups = {'within', 'neighbour', 'random'}; gl = {'same cluster', 'neighbour (0.5-5 s)', 'random (> 60 s)'};
    c = [orange; blue; grey];
    for q = 1:3
        z = S.(groups{q}).rho_z; z = z(~isnan(z));
        if isempty(z); continue; end
        m = mean(z); se = std(z) / sqrt(numel(z));
        bar(ax, q, m, 0.6, 'FaceColor', c(q, :), 'EdgeColor', 'none');
        errorbar(ax, q, m, se, 'k', 'CapSize', 0, 'LineWidth', 1);
        text(ax, q, m + se + 0.03, sprintf('n=%d\ncos %.2f', numel(z), median(S.(groups{q}).cosine, 'omitnan')), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    yline(ax, 0, 'k-');
    set(ax, 'XTick', 1:3, 'XTickLabel', gl);
    ylabel(ax, 'firing-order similarity (z vs neuron shuffle)');
    zw = S.within.rho_z(~isnan(S.within.rho_z)); zr = S.random.rho_z(~isnan(S.random.rho_z));
    pv = NaN; if numel(zw) > 3 && numel(zr) > 3; pv = ranksum(zw, zr); end
    title(ax, sprintf('%s: do ripples share content? (same cluster vs random p = %.2g)', opts.populations{p}, pv), ...
        'FontWeight', 'normal', 'FontSize', 8.5);
    style(ax);
end
end
