function C = classify_ripple_clusters(onset, offset, opts)
% CLASSIFY_RIPPLE_CLUSTERS  Solo / clustered / intermediate ripples after Voroslakos et al.
% (2026, bioRxiv 10.64898/2026.03.27.714843): inter-ripple interval (onset -> onset)
% < 180 ms chains ripples into a cluster (cSPW-R); a ripple >= 500 ms from both
% neighbours is solo (sSPW-R); anything else is intermediate. Decoding windows follow the
% project rule: each ripple [onset, max(offset, onset + 100 ms)]; a cluster from its first
% onset to the end of its last ripple's window.
%
%   C = classify_ripple_clusters(onset, offset)            % defaults below
%   C = classify_ripple_clusters(onset, offset, opts)
%
% opts.cluster_iei (0.18 s), opts.solo_iei (0.5 s), opts.min_window (0.1 s)
% C.type        "solo" | "cluster" | "intermediate" per ripple (sorted by onset)
% C.cluster     cluster id per ripple (0 for non-cluster ripples)
% C.position    position within its cluster (1, 2, ...; 0 otherwise)
% C.window      [n x 2] per-ripple decoding window
% C.clusters    struct: onset, offset (end of last window), members (cell), size
% C.order       sort index applied to the inputs
% C.gmm         two-Gaussian fit to log10 inter-ripple intervals: mu, sigma, weight,
%               intersection (s) between the two means (their boundary was 177 ms)

if nargin < 3; opts = struct(); end
def.cluster_iei = 0.18; def.solo_iei = 0.5; def.min_window = 0.1; def.fit_gmm = true;
for f = fieldnames(def)'
    if ~isfield(opts, f{1}); opts.(f{1}) = def.(f{1}); end
end
[onset, o] = sort(onset(:));
offset = offset(o);
offset = offset(:);
n = numel(onset);
iei = diff(onset);
prev = [inf; iei];
next = [iei; inf];

C.order = o;
C.onset = onset;
C.offset = offset;
C.iei_prev = prev;
C.iei_next = next;
C.type = repmat("intermediate", n, 1);
in_cl = prev < opts.cluster_iei | next < opts.cluster_iei;
C.type(in_cl) = "cluster";
C.type(prev >= opts.solo_iei & next >= opts.solo_iei) = "solo";
C.window = [onset, max(offset, onset + opts.min_window)];

% clusters: maximal chains with consecutive IEI < cluster_iei
C.cluster = zeros(n, 1);
C.position = zeros(n, 1);
k = 0;
i = 1;
members = {};
while i <= n
    if in_cl(i) && next(i) < opts.cluster_iei
        j = i;
        while j < n && next(j) < opts.cluster_iei; j = j + 1; end
        k = k + 1;
        C.cluster(i:j) = k;
        C.position(i:j) = (1:j - i + 1)';
        members{k, 1} = (i:j)'; %#ok<AGROW>
        i = j + 1;
    else
        i = i + 1;
    end
end
C.clusters.members = members;
C.clusters.size = cellfun(@numel, members);
C.clusters.onset = cellfun(@(m) onset(m(1)), members);
C.clusters.offset = cellfun(@(m) max(C.window(m, 2)), members);

C.gmm = struct('mu', NaN(1, 2), 'sigma', NaN(1, 2), 'weight', NaN(1, 2), 'intersection', NaN);
if opts.fit_gmm && numel(iei) > 50
    x = log10(iei(iei > 0));
    try
        g = fitgmdist(x, 2, 'Replicates', 5, 'Options', statset('MaxIter', 500));
        [mu, s] = sort(g.mu);
        sg = sqrt(squeeze(g.Sigma));
        sg = sg(s);
        w = g.ComponentProportion(s);
        xs = linspace(mu(1), mu(2), 2000);
        d = w(1) * normpdf(xs, mu(1), sg(1)) - w(2) * normpdf(xs, mu(2), sg(2));
        ix = find(diff(sign(d)) ~= 0, 1);
        C.gmm = struct('mu', 10.^mu(:)', 'sigma', sg(:)', 'weight', w(:)', 'intersection', 10^xs(max(ix, 1)));
    catch
    end
end
C.opts = opts;
end
