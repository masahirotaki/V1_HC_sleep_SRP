function plot_mean_sem(ax, t, x, col)
% PLOT_MEAN_SEM  Mean (line) +/- SEM (shaded) of the rows of x.
m = mean(x, 1);
s = std(x, 0, 1) / sqrt(size(x, 1));
fill(ax, [t fliplr(t)], [m + s fliplr(m - s)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold(ax, 'on');
plot(ax, t, m, 'Color', col, 'LineWidth', 1.2);
xline(ax, 0, ':');
set(ax, 'TickDir', 'out', 'Box', 'off');
end
