function plot_event_rate_by_state(ax, event_states, states)
% PLOT_EVENT_RATE_BY_STATE  Events per minute spent in each state (count on top of bars).
c = state_colors();
names = ["SWS", "REM", "quietWake", "movement"];
rate = zeros(size(names));
n = zeros(size(names));
for k = 1:numel(names)
    n(k) = sum(event_states == names(k));
    minutes = sum(diff(states.(names(k)), [], 2)) / 60;
    rate(k) = n(k) / max(minutes, eps);
end
b = bar(ax, rate, 'FaceColor', 'flat');
b.CData = cell2mat(arrayfun(@(s) c.(s), names', 'UniformOutput', false));
text(ax, 1:numel(names), rate, cellstr(compose('n=%d', n)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
set(ax, 'XTickLabel', names, 'XTickLabelRotation', 30, 'TickDir', 'out', 'Box', 'off');
ylabel(ax, 'Events / min in state');
end
