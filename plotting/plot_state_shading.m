function plot_state_shading(ax, states)
% PLOT_STATE_SHADING  Shade behavioural states behind the data already in AX (call last).
c = state_colors();
names = ["SWS", "REM", "quietWake", "movement"];
yl = ylim(ax);
hold(ax, 'on');
h = gobjects(1, numel(names));
for k = 1:numel(names)
    iv = states.(names(k));
    h(k) = patch(ax, NaN(1, 4), NaN(1, 4), c.(names(k)), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    for i = 1:size(iv, 1)
        p = patch(ax, iv(i, [1 2 2 1]), yl([1 1 2 2]), c.(names(k)), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        uistack(p, 'bottom');
    end
end
ylim(ax, yl);
xlim(ax, states.recording_range);
legend(ax, h, names, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
end
