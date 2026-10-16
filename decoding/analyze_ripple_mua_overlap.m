function results = analyze_ripple_mua_overlap(subject, date, rec_name, cfg, opts)
% ANALYZE_RIPPLE_MUA_OVERLAP  Do CA1 MUA bursts (1 ms bins, Gaussian sigma 5 ms, z > 3 peak,
% edges at z = 1, <= 300 ms; detect_mua_bursts) and ripples (FindRipples_masa) pick out
% the same events? Also classifies ripples as solo / cluster / intermediate after
% Voroslakos et al. 2026 (classify_ripple_clusters) and fits the two-Gaussian model of
% log inter-ripple intervals.
%
%   results = analyze_ripple_mua_overlap('M00069', '20250929', 'Sleep_Box_1', project_config())
%   results = analyze_ripple_mua_overlap('M00069', '20250929', 'GAVNIK_ABCD_1', project_config())
%
% Ripples: <rec>\detection\extracted_ripple_events.mat (run_sleep_detection) if present,
% else <rec>\replay\replay_<rec>_results.mat (run_abcd_awake_replay).
% Sleep recordings (behavioural_state.mat with SWS): events and the MUA z baseline are
% restricted to SWS; otherwise the whole recording is used.
% Output: <rec>\replay\replay_<rec>_ripple_mua_overlap.mat / .fig / .pdf (vector)

if nargin < 4 || isempty(cfg); cfg = project_config(); end
if nargin < 5; opts = struct(); end
def.state = 'auto';     % 'auto' (SWS for sleep recordings), 'all', or a state name
def.save = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
set(groot, 'defaultFigureVisible', cfg.plot.visible);
paths  = get_session_paths(subject, date, rec_name, cfg);
layout = get_probe_layout(paths, cfg);
units  = load_session_clusters(paths, layout, cfg);
tag = sprintf('%s %s %s', subject, date, strrep(rec_name, '_', ' '));

f_det = fullfile(paths.detection_dir, 'extracted_ripple_events.mat');
if isfile(f_det)
    S = load(f_det); rip = S.ripples;
else
    S = load(fullfile(paths.analysis_dir, 'replay', sprintf('replay_%s_results.mat', rec_name))); rip = S.results.ripples;
end
ca1 = units.cluster_id(units.region == "CA1");
sp = units.spike_times(ismember(units.spike_id, ca1));
t_end = max([sp; rip.offset(:)]) + 1;

state = opts.state;
iv = [0 t_end];
f_state = fullfile(paths.output_dir, 'behavioural_state.mat');
if strcmp(state, 'auto')
    state = 'all';
    if startsWith(rec_name, 'Sleep') && isfile(f_state); state = 'SWS'; end
end
if ~strcmp(state, 'all')
    B0 = load(f_state); iv = B0.behavioural_state.(state);
end
inside = @(t) any(t(:) >= iv(:, 1)' & t(:) < iv(:, 2)', 2);

C = classify_ripple_clusters(rip.onset(:), rip.offset(:));
base_iv = [];
if ~strcmp(state, 'all'); base_iv = iv; end
[Bu, mua] = detect_mua_bursts(sp, t_end, struct('baseline', base_iv));
[Bw, ~] = detect_mua_bursts(sp, t_end);                         % whole-recording baseline, for comparison
r_in = inside(C.onset);
b_in = inside(Bu.peak_time);
bw_in = inside(Bw.peak_time);

ov = @(a0, a1, b0, b1) arrayfun(@(i) any(b0 < a1(i) & b1 > a0(i)), (1:numel(a0))');
rip_has_burst = ov(C.window(:, 1), C.window(:, 2), Bu.onset, Bu.offset);
burst_has_rip = ov(Bu.onset, Bu.offset, C.window(:, 1), C.window(:, 2));
rip_has_burst_w = ov(C.window(:, 1), C.window(:, 2), Bw.onset, Bw.offset);
burst_has_rip_w = ov(Bw.onset, Bw.offset, C.window(:, 1), C.window(:, 2));
pkz = arrayfun(@(i) max(mua.z(mua.t >= C.window(i, 1) & mua.t <= C.window(i, 2))), (1:numel(C.onset))');

types = ["solo", "cluster", "intermediate"];
fprintf('\n%s (%s): %d ripples, %d MUA bursts (+%d longer than 300 ms excluded)\n', rec_name, state, sum(r_in), ...
    sum(b_in), Bu.n_too_long);
fprintf('  ripple types: %s\n', strjoin(arrayfun(@(t) sprintf('%s %.0f%%', t, 100 * mean(C.type(r_in) == t)), types, ...
    'UniformOutput', false), ', '));
fprintf('  IEI two-Gaussian intersection: %.0f ms (fast peak %.0f ms, slow peak %.2f s)\n', 1000 * C.gmm.intersection, ...
    1000 * C.gmm.mu(1), C.gmm.mu(2));
fprintf('  ripples overlapping a MUA burst: %.0f%% (solo %.0f%%, cluster %.0f%%, intermediate %.0f%%)\n', ...
    100 * mean(rip_has_burst(r_in)), 100 * mean(rip_has_burst(r_in & C.type == "solo")), ...
    100 * mean(rip_has_burst(r_in & C.type == "cluster")), 100 * mean(rip_has_burst(r_in & C.type == "intermediate")));
fprintf('  MUA bursts overlapping a ripple: %.0f%%\n', 100 * mean(burst_has_rip(b_in)));
fprintf('  (z over the whole recording instead: ripples with burst %.0f%%, bursts with ripple %.0f%%, n bursts %d)\n', ...
    100 * mean(rip_has_burst_w(r_in)), 100 * mean(burst_has_rip_w(bw_in)), sum(bw_in));
fprintf('  peak MUA z inside the ripple window: median %.1f; > 3 in %.0f%%\n', median(pkz(r_in)), 100 * mean(pkz(r_in) > 3));

results = struct('subject', subject, 'date', date, 'recording', rec_name, 'state', state, 'clusters', C, ...
    'bursts', Bu, 'bursts_wholebaseline', Bw, 'ripple_in_state', r_in, 'burst_in_state', b_in, ...
    'ripple_has_burst', rip_has_burst, 'burst_has_ripple', burst_has_rip, 'ripple_peak_mua_z', pkz);
fig = plot_overlap(results, mua, tag);
if ~opts.save; return; end
out_dir = fullfile(paths.analysis_dir, 'replay');
base = sprintf('replay_%s_ripple_mua_overlap', rec_name);
safe_save(fullfile(out_dir, [base '.mat']), struct('results', results), cfg);
save_vector_figures(fig, {base}, out_dir, cfg);
end


function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8, 'Color', 'none');
end

function f = plot_overlap(r, mua, tag)
C = r.clusters; ri = r.ripple_in_state; bi = r.burst_in_state;
blue = [0.165 0.471 0.839]; orange = [0.922 0.408 0.204]; grey = [0.55 0.55 0.55];
f = figure('Color', 'w', 'Position', [40 40 1600 820]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['%s (%s) - ripple clusters (Voroslakos et al. 2026: IEI < 180 ms = cluster, >= 500 ms = solo) ' ...
    'and CA1 MUA bursts (1 ms, \\sigma = 5 ms, peak z > 3, edges z = 1, <= 300 ms)'], tag, r.state), 'FontSize', 10.5);

ax = nexttile(tl); hold(ax, 'on');
pair_in = ri(1:end - 1) & ri(2:end);
x = log10(diff(C.onset)); x = x(pair_in);
e = linspace(-1.6, 2, 60);
histogram(ax, x, e, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
if ~isnan(C.gmm.intersection)
    xs = linspace(e(1), e(end), 400);
    gcol = [orange; blue];
    for k = 1:2
        plot(ax, xs, C.gmm.weight(k) * normpdf(xs, log10(C.gmm.mu(k)), C.gmm.sigma(k)), '-', 'LineWidth', 2, ...
            'Color', gcol(k, :));
    end
    xline(ax, log10(C.gmm.intersection), 'k-', sprintf('%.0f ms', 1000 * C.gmm.intersection), 'FontSize', 7);
end
xline(ax, log10([0.18 0.5]), 'k:', {'180 ms', '500 ms'}, 'FontSize', 7);
set(ax, 'XTick', -1:2, 'XTickLabel', {'0.1', '1', '10', '100'});
xlabel(ax, 'inter-ripple interval, onset \rightarrow onset (s, log)'); ylabel(ax, 'density');
title(ax, sprintf('fast component peak %.0f ms, slow %.2f s', 1000 * C.gmm.mu(1), C.gmm.mu(2)), 'FontWeight', 'normal');
style(ax);

ax = nexttile(tl); hold(ax, 'on');
sz = C.clusters.size; cl_in = false(numel(sz), 1);
for k = 1:numel(sz); cl_in(k) = ri(C.clusters.members{k}(1)); end
vals = [sum(ri & C.type == "solo"), sum(ri & C.type == "intermediate"), sum(sz(cl_in) == 2), sum(sz(cl_in) == 3), sum(sz(cl_in) >= 4)];
b = bar(ax, vals, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = [blue; grey; orange; orange; orange];
text(ax, b.XEndPoints, vals, string(vals), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
set(ax, 'XTick', 1:5, 'XTickLabel', {'solo ripples', 'intermediate', 'doublets', 'triplets', '4+'});
ylabel(ax, 'count');
title(ax, sprintf('%.0f%% solo, %.0f%% in clusters, %.0f%% intermediate (ripples)', 100 * mean(C.type(ri) == "solo"), ...
    100 * mean(C.type(ri) == "cluster"), 100 * mean(C.type(ri) == "intermediate")), 'FontWeight', 'normal');
style(ax);

ax = nexttile(tl); hold(ax, 'on');
pw = -0.3:0.002:0.5;
grp = {ri & C.type == "solo", ri & C.position == 1};
lab = {'solo ripple', 'first ripple of a cluster'}; col = {blue, orange};
for g = 1:2
    t0 = C.onset(grp{g});
    if isempty(t0); continue; end
    M = interp1(mua.t, mua.z, t0 + pw, 'linear', NaN);
    if size(M, 2) ~= numel(pw); M = M'; end
    mu = mean(M, 1, 'omitnan'); se = std(M, 0, 1, 'omitnan') / sqrt(size(M, 1));
    fill(ax, 1000 * [pw fliplr(pw)], [mu - se fliplr(mu + se)], col{g}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, 1000 * pw, mu, '-', 'Color', col{g}, 'LineWidth', 2, 'DisplayName', sprintf('%s (n = %d)', lab{g}, numel(t0)));
end
yline(ax, [1 3], 'k:', {'edge z = 1', 'peak z = 3'}, 'FontSize', 7, 'HandleVisibility', 'off');
xline(ax, 0, 'k-', 'HandleVisibility', 'off');
xlabel(ax, 'time from ripple onset (ms)'); ylabel(ax, 'CA1 MUA z (\sigma = 5 ms)');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, 'CA1 MUA around ripples', 'FontWeight', 'normal');
style(ax);

ax = nexttile(tl); hold(ax, 'on');
v = 100 * [mean(r.ripple_has_burst(ri)), mean(r.ripple_has_burst(ri & C.type == "solo")), ...
    mean(r.ripple_has_burst(ri & C.type == "cluster")), mean(r.ripple_has_burst(ri & C.type == "intermediate")), ...
    mean(r.burst_has_ripple(bi))];
b = bar(ax, v, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = [grey; blue; orange; grey; [0.3 0.3 0.3]];
text(ax, b.XEndPoints, v, compose('%.0f%%', v), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
set(ax, 'XTick', 1:5, 'XTickLabel', {'all ripples', 'solo', 'cluster', 'intermediate', 'MUA bursts with a ripple'});
ylim(ax, [0 110]); ylabel(ax, '% overlapping');
title(ax, sprintf('ripples containing a MUA burst; %d bursts, %.0f%% contain a ripple', sum(bi), v(5)), 'FontWeight', 'normal');
style(ax);

ax = nexttile(tl); hold(ax, 'on');
histogram(ax, r.ripple_peak_mua_z(ri & C.type == "solo"), 0:0.5:20, 'FaceColor', blue, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'DisplayName', 'solo');
histogram(ax, r.ripple_peak_mua_z(ri & C.type == "cluster"), 0:0.5:20, 'FaceColor', orange, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'DisplayName', 'cluster');
xline(ax, 3, 'k--', 'z = 3', 'HandleVisibility', 'off');
xlabel(ax, 'peak CA1 MUA z inside the ripple window'); ylabel(ax, 'ripples');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, sprintf('median %.1f; %.0f%% above 3', median(r.ripple_peak_mua_z(ri)), 100 * mean(r.ripple_peak_mua_z(ri) > 3)), ...
    'FontWeight', 'normal');
style(ax);

ax = nexttile(tl); hold(ax, 'on');
B = r.bursts;
histogram(ax, 1000 * B.duration(bi & r.burst_has_ripple), 0:10:300, 'FaceColor', orange, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'DisplayName', sprintf('with ripple (n = %d)', sum(bi & r.burst_has_ripple)));
histogram(ax, 1000 * B.duration(bi & ~r.burst_has_ripple), 0:10:300, 'FaceColor', grey, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'DisplayName', sprintf('without ripple (n = %d)', sum(bi & ~r.burst_has_ripple)));
xlabel(ax, 'MUA burst duration (ms)'); ylabel(ax, 'bursts');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, sprintf('burst duration (%d bursts > 300 ms excluded)', B.n_too_long), 'FontWeight', 'normal');
style(ax);
end
