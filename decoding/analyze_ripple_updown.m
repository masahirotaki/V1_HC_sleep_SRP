function results = analyze_ripple_updown(subject, date, sleep_name, cfg, opts)
% ANALYZE_RIPPLE_UPDOWN  Where ripples fall relative to cortical UP / DOWN states (V1 L5,
% DetectSlowWaves_masa): around UP onsets and DOWN onsets, and at the normalised position
% (0 = start, 1 = end) inside each UP and each DOWN state.
%
%   results = analyze_ripple_updown('M00069', '20250929', 'Sleep_Box_2', project_config())
%   results = analyze_ripple_updown(..., struct('stimulus', 'GAVNIK_ABCD_1'))   % + replay groups
%
% Only UP states up to opts.max_up (2 s) are used (longer "UP" periods may not be real UP
% states); DOWN onsets are the ends of those UP states, UP onsets their starts. DOWN
% states in the normalised analysis are the ones that follow a used UP state.
% Ripple groups: all, solo, cluster (all members), first ripple of a cluster (R1),
% intermediate (classify_ripple_clusters, Voroslakos et al. 2026). With opts.stimulus and
% run_sleep_replay results: ripples with significant replay (calibrated, either direction)
% in CA1, in V1, and in both.
% Times are ripple peak times. Output: <sleep>\replay\replay_<sleep>_ripple_updown.mat/.fig/.pdf

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.max_up = 2;
def.peri = [-1 1];
def.bin = 0.02;
def.n_phase_bins = 20;
def.stimulus = '';
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths = get_session_paths(subject, date, sleep_name, cfg);
tag = sprintf('%s %s %s', subject, date, strrep(sleep_name, '_', ' '));
R = load(fullfile(paths.detection_dir, 'extracted_ripple_events.mat')); rip = R.ripples;
SO = load(fullfile(paths.detection_dir, 'extracted_slow_oscillation_events.mat')); so = SO.slow_oscillations;
C = classify_ripple_clusters(rip.onset(:), rip.offset(:));
pk = rip.peaktimes(C.order); pk = pk(:);

groups = struct('name', {'all', 'solo', 'cluster (all)', 'cluster R1', 'intermediate'}, ...
    'mask', {true(size(pk)), C.type == "solo", C.type == "cluster", C.position == 1, C.type == "intermediate"});
if ~isempty(opts.stimulus)
    f = fullfile(paths.analysis_dir, 'replay', sprintf('replay_%s_template_%s_results.mat', sleep_name, opts.stimulus));
    if isfile(f)
        RP = load(f); RP = RP.results;
        sig = struct();
        for p = 1:numel(RP.scores)
            sc = RP.scores(p);
            e = RP.events(sc.event, :);
            is_r = e.kind == "ripple";
            m = false(size(pk));
            m(e.ripple_index(is_r & (sc.forward_cal | sc.reverse_cal))) = true;
            sig.(sc.name) = m;
            groups(end + 1) = struct('name', sprintf('%s replay', strrep(sc.name, '_', ' ')), 'mask', m); %#ok<AGROW>
        end
        f2 = fieldnames(sig);
        if numel(f2) >= 2
            groups(end + 1) = struct('name', 'CA1 & V1 replay', 'mask', sig.(f2{1}) & sig.(f2{2}));
        end
    else
        warning('No replay results %s: replay groups skipped.', f);
    end
end

up = so.UP; dn = so.DOWN;
up_dur = diff(up, [], 2);
ok_up = up_dur > 0 & up_dur <= opts.max_up;
U = up(ok_up, :);
% DOWN states that start where a used UP state ends
[~, idn] = min(abs(dn(:, 1) - U(:, 2)'), [], 1);
follow = abs(dn(idn, 1) - U(:, 2)) < 0.01;
D = dn(idn(follow), :);
fprintf('%s: %d UP states (%d <= %.1f s used), %d DOWN states (%d following a used UP); %d ripples\n', sleep_name, ...
    size(up, 1), size(U, 1), opts.max_up, size(dn, 1), size(D, 1), numel(pk));

edges = opts.peri(1):opts.bin:opts.peri(2);
ctr = edges(1:end - 1) + opts.bin / 2;
pedges = linspace(0, 1, opts.n_phase_bins + 1);
for g = 1:numel(groups)
    t = pk(groups(g).mask);
    groups(g).n = numel(t);
    groups(g).peri_up = peri(t, U(:, 1));          % ripples / s around UP onsets
    groups(g).peri_down = peri(t, U(:, 2));        % ... around DOWN onsets
    [groups(g).phase_up, groups(g).n_in_up] = phase(t, U);
    [groups(g).phase_down, groups(g).n_in_down] = phase(t, D);
end
results = struct('subject', subject, 'date', date, 'sleep', sleep_name, 'opts', opts, 'UP', U, 'DOWN', D, ...
    'UP_all', up, 'DOWN_all', dn, 'peri_t', ctr, 'phase_edges', pedges, 'groups', groups, 'clusters', C);
for g = 1:numel(groups)
    fprintf('  %-16s n = %4d | in UP %3d (%.0f%% of group), in DOWN %3d\n', groups(g).name, groups(g).n, ...
        groups(g).n_in_up, 100 * groups(g).n_in_up / max(groups(g).n, 1), groups(g).n_in_down);
end
fig = plot_updown(results, tag);
if ~opts.save; return; end
out_dir = fullfile(paths.analysis_dir, 'replay');
base = sprintf('replay_%s_ripple_updown', sleep_name);
safe_save(fullfile(out_dir, [base '.mat']), struct('results', results), cfg);
save_vector_figures(fig, {base}, out_dir, cfg);

    function r = peri(t, t0)
        % mean rate (events / s) of ripple peaks around each reference time
        c = zeros(1, numel(edges) - 1);
        for k = 1:numel(t0)
            d = t - t0(k);
            d = d(d >= edges(1) & d < edges(end));
            if ~isempty(d); c = c + histcounts(d, edges); end
        end
        r = c / (numel(t0) * opts.bin);
    end

    function [h, n_in] = phase(t, iv)
        % ripples per state per unit normalised time, and how many ripples fall inside
        ph = [];
        for k = 1:size(iv, 1)
            in = t >= iv(k, 1) & t < iv(k, 2);
            ph = [ph; (t(in) - iv(k, 1)) / (iv(k, 2) - iv(k, 1))]; %#ok<AGROW>
        end
        n_in = numel(ph);
        h = histcounts(ph, pedges) / (size(iv, 1) * diff(pedges(1:2)));
    end
end


function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function f = plot_updown(r, tag)
G = r.groups;
cols = [0.15 0.15 0.15; 0.165 0.471 0.839; 0.922 0.408 0.204; 0.65 0.25 0.1; 0.55 0.55 0.55; ...
    0.106 0.686 0.478; 0.290 0.227 0.655; 0.8 0.6 0.0];
show_peri = 1:min(numel(G), 8);
f = figure('Color', 'w', 'Position', [40 40 1600 860]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s - ripples (peak times) around V1 UP / DOWN states; UP states <= %.1f s (%d UP, %d DOWN); ' ...
    'solo / cluster after Voroslakos et al. 2026'], tag, r.opts.max_up, size(r.UP, 1), size(r.DOWN, 1)), 'FontSize', 10.5);
t = 1000 * r.peri_t;
ttl = {'UP onset (DOWN \rightarrow UP)', 'DOWN onset (UP \rightarrow DOWN)'};
fld = {'peri_up', 'peri_down'};
for q = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for g = show_peri
        y = G(g).(fld{q});
        if G(g).n == 0; continue; end
        plot(ax, t, movmean(y, 3), '-', 'Color', cols(g, :), 'LineWidth', 1.6 + 0.6 * (g == 1), ...
            'DisplayName', sprintf('%s (n = %d)', G(g).name, G(g).n));
    end
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    xlabel(ax, sprintf('time from %s (ms)', ttl{q})); ylabel(ax, 'ripple rate (Hz)');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
    title(ax, ['ripples around ' ttl{q}], 'FontWeight', 'normal');
    style(ax);
end
ax = nexttile(tl); hold(ax, 'on');
histogram(ax, diff(r.UP_all, [], 2), 0:0.1:6, 'FaceColor', [0.165 0.471 0.839], 'EdgeColor', 'none', 'DisplayName', 'UP');
histogram(ax, diff(r.DOWN_all, [], 2), 0:0.1:6, 'FaceColor', [0.922 0.408 0.204], 'EdgeColor', 'none', 'DisplayName', 'DOWN');
xline(ax, r.opts.max_up, 'k--', sprintf('UP cut %.0f s', r.opts.max_up), 'FontSize', 7, 'HandleVisibility', 'off');
xlabel(ax, 'duration (s)'); ylabel(ax, 'states');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, sprintf('%.0f%% of UP states are <= %.0f s', 100 * mean(diff(r.UP_all, [], 2) <= r.opts.max_up), r.opts.max_up), ...
    'FontWeight', 'normal');
style(ax);
pc = (r.phase_edges(1:end - 1) + r.phase_edges(2:end)) / 2;
pf = {'phase_up', 'phase_down'}; pn = {'n_in_up', 'n_in_down'}; pt = {'inside UP states', 'inside DOWN states'};
for q = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for g = show_peri
        if G(g).n == 0; continue; end
        y = G(g).(pf{q});
        plot(ax, pc, y, '-o', 'Color', cols(g, :), 'MarkerSize', 3, 'MarkerFaceColor', cols(g, :), ...
            'LineWidth', 1.4 + 0.6 * (g == 1), 'DisplayName', sprintf('%s (%d inside)', G(g).name, G(g).(pn{q})));
    end
    xlabel(ax, sprintf('normalised position %s (0 = start, 1 = end)', pt{q})); ylabel(ax, 'ripples per state per unit time');
    legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
    title(ax, ['ripple position ' pt{q}], 'FontWeight', 'normal');
    style(ax);
end
ax = nexttile(tl); hold(ax, 'on');
ud = diff(r.UP, [], 2);
nin = arrayfun(@(k) sum(r.clusters.window(:, 1) >= r.UP(k, 1) & r.clusters.window(:, 1) < r.UP(k, 2)), 1:size(r.UP, 1));
eb = 0:0.25:r.opts.max_up;
mr = arrayfun(@(i) mean(nin(ud >= eb(i) & ud < eb(i + 1))), 1:numel(eb) - 1);
bar(ax, eb(1:end - 1) + 0.125, mr, 1, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'w');
xlabel(ax, 'UP state duration (s)'); ylabel(ax, 'ripples per UP state');
title(ax, 'ripples per UP state vs its duration', 'FontWeight', 'normal');
style(ax);
end
