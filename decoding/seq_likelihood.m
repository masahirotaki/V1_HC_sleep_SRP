function R = seq_likelihood(t, id, W, n_shuffle, rng_seed, do_contrib)
% SEQ_LIKELIHOOD  Log-sequence likelihood Lseq of one spike sequence under an
% ordered firing probability matrix (seq_order_matrix), with a cell-identity shuffle
% null, as in Huh, Yun, ..., Jung (2026) Nat Commun (Estimate_replay_pMat_f1.m).
%
%   R = seq_likelihood(t, id, W, n_shuffle, rng_seed)
%
% t, id   spike times and neuron indices of the event / window (any order)
% W       seq_order_matrix(...).W = log p - log p'  (W(i, j) > 0: i tends to fire after j)
% Lseq = sum over consecutive spikes W(later, earlier) + sum over all spike pairs
%        W(later, earlier); pairs at identical times are skipped. Same-neuron pairs
%        contribute 0. Positive = forward (template order), negative = reverse.
%        (= their x_p_replay; the paper's eq. 12 differs by a factor 1/2.)
% Null: neuron identities permuted across the spikes (times kept), n_shuffle times.
%
% R.L, R.L_shuffle, R.n_pairs (contributing pairs), R.L_norm = L / n_pairs,
% R.z = (L - mean(null)) / sd(null), R.pct = fraction of null below L (ties count 1/2),
% R.n_spikes, R.n_neurons
%
% do_contrib (default false): per-neuron decomposition. Each pair term is split half to
% each of its two neurons (consecutive-pair terms likewise), so sum(R.contrib) = R.L.
% R.contrib       [n x 1] neuron contributions to L (0 for neurons without spikes)
% R.contrib_null  [n x 1] mean contribution over the shuffles
% R.contrib_z     [n x 1] (contrib - contrib_null) / sd(null): sums to R.z

if nargin < 4; n_shuffle = 1000; end
if nargin >= 5 && ~isempty(rng_seed); rng(rng_seed); end
if nargin < 6; do_contrib = false; end
[t, o] = sort(t(:));
id = id(o);
id = id(:);
n = size(W, 1);
N = numel(t);
R.n_spikes = N;
R.n_neurons = numel(unique(id));
if N < 2
    [R.L, R.n_pairs, R.L_norm, R.z, R.pct] = deal(NaN);
    R.L_shuffle = [];
    if do_contrib; [R.contrib, R.contrib_null, R.contrib_z] = deal(zeros(n, 1)); end
    return
end
[~, ~, g] = unique(t);
Gt = sparse(g, 1:N, 1, max(g), N);
diff_time = diff(g) ~= 0;

if do_contrib
    [R.L, R.n_pairs, R.contrib] = core(id);
else
    [R.L, R.n_pairs] = core(id);
end
R.L_norm = R.L / R.n_pairs;
R.L_shuffle = zeros(n_shuffle, 1);
c_null = zeros(n, 1);
for k = 1:n_shuffle
    if do_contrib
        [R.L_shuffle(k), ~, c] = core(id(randperm(N)));
        c_null = c_null + c;
    else
        R.L_shuffle(k) = core(id(randperm(N)));
    end
end
sd = std(R.L_shuffle);
R.z = (R.L - mean(R.L_shuffle)) / sd;
R.pct = mean(R.L_shuffle < R.L) + 0.5 * mean(R.L_shuffle == R.L);
if do_contrib
    R.contrib_null = c_null / n_shuffle;
    R.contrib_z = (R.contrib - R.contrib_null) / sd;
end

    function [L, np, c] = core(ids)
        X = sparse(1:N, ids, 1, N, n);
        G = Gt * X;
        prev = [zeros(1, n); cumsum(full(G(1:end - 1, :)), 1)];
        Hprev = prev(g, :);                                  % earlier spikes per neuron
        PW = W(ids, :) .* Hprev;                             % row: later spike, column: earlier neuron
        cons = W(sub2ind([n n], ids(2:end), ids(1:end - 1))) .* diff_time;   % consecutive pairs
        L = sum(PW(:)) + sum(cons);
        np = [];
        if nargout > 1
            np = sum(sum(Hprev, 2) - Hprev(sub2ind([N n], (1:N)', ids)));
        end
        if nargout > 2   % half of every pair term to each of its two neurons
            c = 0.5 * (accumarray(ids, sum(PW, 2), [n 1]) + sum(PW, 1)') ...
                + 0.5 * (accumarray(ids(2:end), cons, [n 1]) + accumarray(ids(1:end - 1), cons, [n 1]));
        end
    end
end
