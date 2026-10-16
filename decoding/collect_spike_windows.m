function [w, n_sp] = collect_spike_windows(spike_times, spike_id, ids, starts, stops)
% COLLECT_SPIKE_WINDOWS  Spikes of the units IDS inside each window [starts(k), stops(k)).
%
% w{k}  [time re starts(k), neuron index into ids] (sorted by time)
% n_sp  spikes per unit summed over all windows

in = ismember(spike_id, ids);
st = spike_times(in);
[~, col] = ismember(spike_id(in), ids);
[st, o] = sort(st);
col = col(o);
w = cell(numel(starts), 1);
n_sp = zeros(numel(ids), 1);
for k = 1:numel(starts)
    a = find(st >= starts(k), 1);
    b = find(st < stops(k), 1, 'last');
    if isempty(a) || isempty(b) || b < a
        w{k} = zeros(0, 2);
        continue
    end
    w{k} = [st(a:b) - starts(k), col(a:b)];
    n_sp = n_sp + accumarray(col(a:b), 1, [numel(ids) 1]);
end
end
