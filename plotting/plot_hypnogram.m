function plot_hypnogram(ax, states)
% PLOT_HYPNOGRAM  State intervals as coloured bars: movement (top) to REM (bottom).
c = state_colors();
names = ["REM", "SWS", "quietWake", "movement"];
hold(ax, 'on');
for k = 1:numel(names)
    iv = states.(names(k));
    for i = 1:size(iv, 1)
        patch(ax, iv(i, [1 2 2 1]), k + [-0.4 -0.4 0.4 0.4], c.(names(k)), 'EdgeColor', 'none');
    end
end
set(ax, 'YTick', 1:numel(names), 'YTickLabel', names, 'YLim', [0.4 numel(names) + 0.6], ...
    'TickDir', 'out', 'Box', 'off');
xlim(ax, states.recording_range);
end
