function fig = plot_sleep_hypnogram(states, cfg, title_str)
% PLOT_SLEEP_HYPNOGRAM  Sleep/wake states over time, one row per state
% (Movement, Quiet Wake, SWS, REM), drawn as in Sleep_analysis_Ellie.

rows = {'movement', 'Movement', 'k', 4; 'quietWake', 'Quiet Wake', 'g', 3; ...
        'SWS', 'SWS', 'b', 2; 'REM', 'REM', 'r', 1};
fig = figure('Name', 'Sleep stages', 'Position', [100 100 1400 500], 'Visible', cfg.plot.visible);
ax = axes(fig, 'Position', [0.12 0.14 0.85 0.76]); % left margin for the state labels
hold(ax, 'on');
for r = 1:size(rows, 1)
    iv = states.(rows{r, 1});
    y = rows{r, 4} + [-0.5 0.5];
    for i = 1:size(iv, 1)
        patch(ax, iv(i, [1 2 2 1]), y([1 1 2 2]), rows{r, 3}, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end
end
xlim(ax, states.recording_range);
ylim(ax, [0 5]);
yticks(ax, 1:4);
yticklabels(ax, {});
for r = 1:size(rows, 1)
    text(ax, ax.XLim(1) - 0.01 * diff(ax.XLim), rows{r, 4}, rows{r, 2}, 'Color', rows{r, 3}, ...
        'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
xlabel(ax, 'Time (s)');
set(ax, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14);
title(ax, sprintf('%s | Sleep/Wake States over time', title_str), 'Interpreter', 'none', 'FontSize', 16);
end
