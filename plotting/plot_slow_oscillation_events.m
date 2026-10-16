function fig = plot_slow_oscillation_events(so, l5_lfp, tvec, fs, units, states, cfg, title_str)
% PLOT_SLOW_OSCILLATION_EVENTS  UP/DOWN detection summary: SO rate over time, DOWN-peak
% triggered L5 LFP and V1 MUA (DOWN states should be silent), UP/DOWN durations and
% an example SWS segment with the MUA raster.

fig = figure('Name', 'Slow oscillations', 'Position', [50 50 1500 1000], 'Visible', cfg.plot.visible);
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact');
title(tl, sprintf('%s | UP/DOWN on L5, shank %d depth %g um | %d DOWN states | %d V1 units', ...
    title_str, so.channel.shank, so.channel.depth, size(so.DOWN, 1), numel(so.mua_cluster_id)), 'Interpreter', 'none');
if isempty(so.peaktimes)
    text(nexttile(tl), 0.5, 0.5, 'No slow oscillations detected', 'HorizontalAlignment', 'center');
    return
end
delta = bandpass_lfp(l5_lfp, fs, cfg.so.filterparms.deltafilter);
in_mua = ismember(units.spike_id, so.mua_cluster_id);
mua_times = units.spike_times(in_mua);
mua_id = units.spike_id(in_mua);

% SO rate over the session
ax = nexttile(tl, [1 2]);
edges = states.recording_range(1):30:states.recording_range(2);
counts = histcounts(so.peaktimes, edges);
stairs(ax, edges(1:end-1), counts / 0.5, 'k', 'HandleVisibility', 'off');
plot_state_shading(ax, states);
ylabel(ax, 'SO / min (30 s bins)');
xlabel(ax, 'Time (s)');
title(ax, 'Slow-oscillation rate (detected within SWS only; shading = state)');
set(ax, 'TickDir', 'out', 'Box', 'off');

ax = nexttile(tl);
histogram(ax, so.peak_magnitude, 40, 'FaceColor', [0.4 0.4 0.4]);
xlabel(ax, 'SW peak magnitude (modZ)'); ylabel(ax, 'Count'); set(ax, 'TickDir', 'out', 'Box', 'off');

% DOWN-peak triggered L5 LFP and V1 MUA
win = [-1 1];
[raw, t_rel] = event_snippets(l5_lfp, tvec, fs, so.peaktimes, win);
flt = event_snippets(delta, tvec, fs, so.peaktimes, win);
ax = nexttile(tl);
plot_mean_sem(ax, t_rel, raw, [0 0 0]);
hold(ax, 'on');
plot(ax, t_rel, mean(flt, 1), 'Color', [0.2 0.4 0.8], 'LineWidth', 1.2);
xlabel(ax, 'Time from SW (DOWN) peak (s)'); ylabel(ax, '\muV');
title(ax, sprintf('L5 LFP: raw (black), %g-%g Hz (blue)', cfg.so.filterparms.deltafilter));

ax = nexttile(tl);
bin = 0.01;
mua_edges = tvec(1):bin:tvec(end);
mua_rate = histcounts(mua_times, mua_edges) / bin / max(numel(so.mua_cluster_id), 1);
[mua_snips, t_mua] = event_snippets(mua_rate, mua_edges(1:end-1) + bin / 2, 1 / bin, so.peaktimes, win);
plot_mean_sem(ax, t_mua, smoothdata(mua_snips, 2, 'gaussian', 5), [0.85 0.3 0.3]);
xlabel(ax, 'Time from SW (DOWN) peak (s)'); ylabel(ax, 'Spikes / s / unit');
title(ax, 'V1 MUA around DOWN states');

ax = nexttile(tl);
histogram(ax, so.DOWN_duration * 1000, 40, 'FaceColor', [0.2 0.4 0.8]);
xlabel(ax, 'DOWN duration (ms)'); ylabel(ax, 'Count'); set(ax, 'TickDir', 'out', 'Box', 'off');
ax = nexttile(tl);
histogram(ax, log10(so.UP_duration), 40, 'FaceColor', [0.85 0.3 0.3]);
xlabel(ax, 'log_{10} UP duration (s)'); ylabel(ax, 'Count'); set(ax, 'TickDir', 'out', 'Box', 'off');

% Example 10 s of SWS with the densest SO activity
[~, densest] = max(counts);
t0 = max(edges(densest), tvec(1));
seg = [t0 t0 + 10];
k = tvec >= seg(1) & tvec <= seg(2);
ax = nexttile(tl, [1 2]);
hold(ax, 'on');
lim = [min(l5_lfp(k)) max(l5_lfp(k))];
down = so.DOWN(so.DOWN(:, 2) >= seg(1) & so.DOWN(:, 1) <= seg(2), :);
for i = 1:size(down, 1)
    patch(ax, down(i, [1 2 2 1]), lim([1 1 2 2]), [0.2 0.4 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
plot(ax, tvec(k), l5_lfp(k), 'Color', [0.5 0.5 0.5]);
plot(ax, tvec(k), delta(k), 'Color', [0.2 0.4 0.8], 'LineWidth', 1.2);
ylabel(ax, 'L5 LFP (\muV)');
yyaxis(ax, 'right');
in_seg = mua_times >= seg(1) & mua_times <= seg(2);
[~, by_depth] = sort(units.depth);
[~, depth_rank] = ismember(mua_id(in_seg), units.cluster_id(by_depth));
plot(ax, mua_times(in_seg), depth_rank, 'k.', 'MarkerSize', 3);
ylabel(ax, 'V1 unit (by depth)');
xlim(ax, seg);
xlabel(ax, 'Time (s)');
title(ax, 'Example: DOWN states (blue shading), L5 LFP and V1 spikes');
end
