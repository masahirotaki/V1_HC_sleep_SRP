function counts = bin_trial_spikes(spike_times, spike_id, cluster_ids, align, edges)
% BIN_TRIAL_SPIKES  Spike counts around event times.
%
%   counts = bin_trial_spikes(spike_times, spike_id, cluster_ids, align, edges)
%
% align  [nTrials x 1] event times (s); edges  bin edges relative to the event (s)
% counts [nTrials x nNeurons x nBins], neurons in the order of cluster_ids.
% Windows of different trials may overlap.

align = align(:);
n_trials = numel(align);
n_bins = numel(edges) - 1;
counts = zeros(n_trials, numel(cluster_ids), n_bins);

keep = spike_times >= min(align) + edges(1) & spike_times < max(align) + edges(end);
[is_unit, col] = ismember(spike_id(keep), cluster_ids);
st = spike_times(keep);
st = st(is_unit);
col = col(is_unit);
for n = 1:numel(cluster_ids)
    s = st(col == n);
    if isempty(s); continue; end
    rel = s - align';                         % spikes x trials
    in_win = rel >= edges(1) & rel < edges(end);
    [~, trial] = find(in_win);
    bin = discretize(rel(in_win), edges);
    ok = ~isnan(bin(:));
    trial = trial(:); bin = bin(:);              % column vectors even when a neuron has one spike
    counts(:, n, :) = reshape(accumarray([trial(ok) bin(ok)], 1, [n_trials n_bins]), n_trials, 1, n_bins);
end
end
