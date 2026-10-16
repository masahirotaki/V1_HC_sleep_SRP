function fig = plot_oscillation_events(events, lfp, tvec, fs, states, opt, cfg, title_str)
% PLOT_OSCILLATION_EVENTS  Detection summary for ripples or spindles: timeline, rate per
% state, SWS-event-triggered averages, example events and peak / duration /
% inter-event-interval distributions.
%
% opt.name            'Ripple' | 'Spindle'
% opt.band            band-pass used for detection (Hz)
% opt.avg_win         window for triggered averages (s)
% opt.example_win     window for example traces (s)
% opt.peak_threshold  detection peak threshold (z), drawn on the histogram

fig = figure('Name', opt.name, 'Position', [50 50 1500 950], 'Visible', cfg.plot.visible);
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact');
title(tl, sprintf('%s | %ss on shank %d, depth %g um | %d events', title_str, opt.name, ...
    events.channel.shank, events.channel.depth, numel(events.peaktimes)), 'Interpreter', 'none');
if isempty(events.peaktimes)
    text(nexttile(tl), 0.5, 0.5, sprintf('No %ss detected', lower(opt.name)), 'HorizontalAlignment', 'center');
    return
end
c = state_colors();
filtered = bandpass_lfp(lfp, fs, opt.band);
ms = 1000;

% Timeline
ax = nexttile(tl, [1 2]);
hold(ax, 'on');
for s = unique(events.state)'
    k = events.state == s;
    scatter(ax, events.peaktimes(k), events.peak_zscore(k), 8, c.(s), 'filled', 'HandleVisibility', 'off');
end
plot_state_shading(ax, states);
ylabel(ax, 'Peak z');
xlabel(ax, 'Time (s)');
title(ax, sprintf('%s peak z-score over the session (shading = state)', opt.name));
set(ax, 'TickDir', 'out', 'Box', 'off');

ax = nexttile(tl);
plot_event_rate_by_state(ax, events.state, states);
title(ax, sprintf('%s rate per state', opt.name));

% SWS-event triggered averages (all events if none in SWS)
use = events.state == "SWS";
label = 'SWS';
if ~any(use); use = true(size(events.state)); label = 'all'; end
[raw, t_rel] = event_snippets(lfp, tvec, fs, events.peaktimes(use), opt.avg_win);
flt = event_snippets(filtered, tvec, fs, events.peaktimes(use), opt.avg_win);
ax = nexttile(tl);
plot_mean_sem(ax, t_rel * ms, raw, [0 0 0]);
title(ax, sprintf('Mean raw LFP, %s %ss (n=%d)', label, lower(opt.name), size(raw, 1)));
xlabel(ax, 'Time from peak (ms)'); ylabel(ax, '\muV');
ax = nexttile(tl);
plot_mean_sem(ax, t_rel * ms, flt, [0.8 0.2 0.2]);
title(ax, sprintf('Mean %g-%g Hz filtered LFP', opt.band));
xlabel(ax, 'Time from peak (ms)'); ylabel(ax, '\muV');

% Example events spread over the peak-z distribution
ax = nexttile(tl);
[~, order] = sort(events.peak_zscore);
examples = order(round(linspace(1, numel(order), min(8, numel(order)))));
[ex_raw, t_ex] = event_snippets(lfp, tvec, fs, events.peaktimes(examples), opt.example_win);
ex_flt = event_snippets(filtered, tvec, fs, events.peaktimes(examples), opt.example_win);
gap = 1.5 * median(max(ex_raw, [], 2) - min(ex_raw, [], 2));
hold(ax, 'on');
for i = 1:size(ex_raw, 1)
    plot(ax, t_ex * ms, ex_raw(i, :) - mean(ex_raw(i, :)) + i * gap, 'Color', [0.5 0.5 0.5]);
    plot(ax, t_ex * ms, ex_flt(i, :) + i * gap, 'Color', [0.8 0.2 0.2]);
end
set(ax, 'YTick', [], 'TickDir', 'out', 'Box', 'off');
xlabel(ax, 'Time from peak (ms)');
title(ax, 'Examples, raw (grey) and filtered (red), low \rightarrow high peak z');

ax = nexttile(tl);
histogram(ax, events.peak_zscore, 40, 'FaceColor', [0.4 0.4 0.4]);
xline(ax, opt.peak_threshold, 'r--');
xlabel(ax, 'Peak z-score'); ylabel(ax, 'Count'); set(ax, 'TickDir', 'out', 'Box', 'off');
ax = nexttile(tl);
histogram(ax, events.duration * ms, 40, 'FaceColor', [0.4 0.4 0.4]);
xlabel(ax, 'Duration (ms)'); ylabel(ax, 'Count'); set(ax, 'TickDir', 'out', 'Box', 'off');
ax = nexttile(tl);
histogram(ax, log10(diff(events.peaktimes)), 50, 'FaceColor', [0.4 0.4 0.4]);
xlabel(ax, sprintf('log_{10} inter-%s interval (s)', lower(opt.name))); ylabel(ax, 'Count');
set(ax, 'TickDir', 'out', 'Box', 'off');
end
