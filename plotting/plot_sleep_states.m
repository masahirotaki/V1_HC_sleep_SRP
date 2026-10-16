function fig = plot_sleep_states(states, Behaviour, l5_lfp, ca1_lfp, tvec, fs, cfg, title_str)
% PLOT_SLEEP_STATES  Mobility, hypnogram, L5 / CA1 spectrograms and the band-power
% features the sleep classifier relies on (cortical delta, hippocampal theta/delta).

fig = figure('Name', 'Sleep states', 'Position', [50 50 1400 950], 'Visible', cfg.plot.visible);
tl = tiledlayout(fig, 6, 1, 'TileSpacing', 'compact');
title(tl, title_str, 'Interpreter', 'none');

ax(1) = nexttile(tl);
plot(ax(1), Behaviour.tvec, Behaviour.mobility_change, 'k');
hold(ax(1), 'on');
yline(ax(1), cfg.behaviour.change_threshold, 'r--', 'threshold');
moving = double(Behaviour.moving);
moving(~Behaviour.moving) = NaN;
plot(ax(1), Behaviour.tvec, moving * cfg.behaviour.change_threshold * 3, 'r', 'LineWidth', 3);
ylabel(ax(1), '|\Delta mobility|');
set(ax(1), 'YLim', [0 cfg.behaviour.change_threshold * 4], 'TickDir', 'out', 'Box', 'off');
title(ax(1), 'Video mobility change (red bar = moving)');

ax(2) = nexttile(tl);
plot_hypnogram(ax(2), states);
title(ax(2), 'Sleep / wake states');

[l5_spec, f, ts] = lfp_spectrogram(l5_lfp, fs, tvec(1));
[ca1_spec, ~, ~] = lfp_spectrogram(ca1_lfp, fs, tvec(1));
ax(3) = nexttile(tl);
plot_spectrogram(ax(3), ts, f, l5_spec);
title(ax(3), 'L5 LFP spectrogram (dB)');
ax(4) = nexttile(tl);
plot_spectrogram(ax(4), ts, f, ca1_spec);
title(ax(4), 'CA1 LFP spectrogram (dB)');

band = @(spec, lim) mean(spec(f >= lim(1) & f <= lim(2), :), 1);
ax(5) = nexttile(tl);
plot(ax(5), ts, smoothdata(zscore(log(band(l5_spec, [0.5 4]))), 'movmean', 10), 'Color', [0.2 0.4 0.8]);
ylabel(ax(5), 'z'); title(ax(5), 'L5 delta (0.5-4 Hz) power, 10 s smoothed');
set(ax(5), 'TickDir', 'out', 'Box', 'off');
ax(6) = nexttile(tl);
plot(ax(6), ts, smoothdata(band(ca1_spec, [6 10]) ./ band(ca1_spec, [1 4]), 'movmean', 10), 'Color', [0.85 0.3 0.3]);
ylabel(ax(6), 'ratio'); title(ax(6), 'CA1 theta (6-10 Hz) / delta (1-4 Hz), 10 s smoothed');
xlabel(ax(6), 'Time (s)');
set(ax(6), 'TickDir', 'out', 'Box', 'off');

linkaxes(ax, 'x');
xlim(ax(1), states.recording_range);
end


function [p, f, ts] = lfp_spectrogram(x, fs, t0)
% 4 s windows, 1 s steps, 0.5-30 Hz, on a 125 Hz copy of the signal.
r = round(fs / 125);
xd = decimate(double(x(:)), r);
fsd = fs / r;
f = 0.5:0.25:30;
[s, f, ts] = spectrogram(xd, hann(round(4 * fsd)), round(3 * fsd), f, fsd);
p = abs(s) .^ 2;
ts = ts + t0;
end


function plot_spectrogram(ax, ts, f, p)
pdb = 10 * log10(p);
imagesc(ax, ts, f, pdb);
axis(ax, 'xy');
clim(ax, prctile(pdb(:), [5 99]));
colormap(ax, 'parula');
ylabel(ax, 'Hz');
set(ax, 'TickDir', 'out', 'Box', 'off');
end
