function label = burst_labels(st, isi_thr)
% BURST_LABELS  Label each spike as single (0), burst onset (1) or later burst spike (2).
%
% Same rule as Manuscript/bursts.py (detect_bursts, min_spikes = 2): a burst is a
% run of spikes whose consecutive ISIs are all <= isi_thr (s). st must be sorted.
% Fraction of spikes in bursts (FSB) = mean(label > 0).

st = st(:);
label = zeros(size(st));
if numel(st) < 2; return; end
b = diff(st) <= isi_thr & diff(st) > 0;  % interval i joins spikes i and i+1
in_burst = [b; false] | [false; b];
starts = in_burst & ~[false; b];
label(in_burst) = 2;
label(starts) = 1;
end
