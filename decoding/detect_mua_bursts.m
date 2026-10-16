function [B, mua] = detect_mua_bursts(spike_times, t_end, opts)
% DETECT_MUA_BURSTS  Population bursts from pooled spikes (e.g. all CA1 units).
% MUA in 1 ms bins, smoothed with a Gaussian kernel (sigma 5 ms), z-scored; a burst is a
% stretch with z > opts.edge_z (onset/offset at the z = 1 crossings) whose peak z exceeds
% opts.peak_z and whose duration is at most opts.max_duration.
%
%   [B, mua] = detect_mua_bursts(spike_times, t_end)
%   [B, mua] = detect_mua_bursts(spike_times, t_end, struct('baseline', SWS_intervals))
%
% opts.bin (0.001 s), opts.sigma (0.005 s), opts.peak_z (3), opts.edge_z (1),
% opts.max_duration (0.3 s), opts.baseline ([] = whole recording; else [n x 2] intervals
% used for the z-score mean / SD)
% B.onset, B.offset, B.peak_time, B.peak_z, B.duration (column vectors)
% mua.t (bin centres), mua.z

if nargin < 3; opts = struct(); end
def.bin = 0.001; def.sigma = 0.005; def.peak_z = 3; def.edge_z = 1; def.max_duration = 0.3; def.baseline = [];
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
edges = 0:opts.bin:t_end + opts.bin;
t = edges(1:end - 1) + opts.bin / 2;
cnt = histcounts(spike_times, edges);
s = opts.sigma / opts.bin;
k = exp(-0.5 * ((-ceil(4 * s):ceil(4 * s)) / s).^2);
k = k / sum(k);
m = conv(cnt, k, 'same');
use = true(size(t));
if ~isempty(opts.baseline)
    use = false(size(t));
    for i = 1:size(opts.baseline, 1); use(t >= opts.baseline(i, 1) & t < opts.baseline(i, 2)) = true; end
end
z = (m - mean(m(use))) / std(m(use));

above = z > opts.edge_z;
d = diff([false above false]);
st = find(d == 1);
en = find(d == -1) - 1;
B = struct('onset', [], 'offset', [], 'peak_time', [], 'peak_z', [], 'duration', []);
for i = 1:numel(st)
    [pz, j] = max(z(st(i):en(i)));
    dur = (en(i) - st(i) + 1) * opts.bin;
    if pz > opts.peak_z && dur <= opts.max_duration
        B.onset(end + 1, 1) = t(st(i)) - opts.bin / 2;
        B.offset(end + 1, 1) = t(en(i)) + opts.bin / 2;
        B.peak_time(end + 1, 1) = t(st(i) + j - 1);
        B.peak_z(end + 1, 1) = pz;
        B.duration(end + 1, 1) = dur;
    end
end
B.n_too_long = sum(arrayfun(@(i) max(z(st(i):en(i))) > opts.peak_z && (en(i) - st(i) + 1) * opts.bin > opts.max_duration, ...
    1:numel(st)));
B.opts = opts;
mua.t = t;
mua.z = z;
end
