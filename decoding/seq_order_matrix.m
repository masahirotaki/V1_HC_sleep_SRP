function M = seq_order_matrix(windows, n)
% SEQ_ORDER_MATRIX  Ordered firing probability matrix of Huh, Yun, ..., Jung (2026)
% Nat Commun 17:8304 ("A likelihood-based method for identifying replay from spike
% sequences"), ported from their Make_pMat_f1.m (Figshare 10.6084/m9.figshare.31942599).
%
%   M = seq_order_matrix(windows, n)
%
% windows  cell array, one entry per active-phase window (trial): [t id] with spike
%          times (s) and neuron indices 1..n
% M.n_after(i, j)  number of spike pairs (within a window) where neuron i fires AFTER j
% M.n_pairs(i, j)  number of comparable pairs (spikes at identical times excluded)
% M.p_raw          n_after ./ n_pairs, clamped to [0.01 0.99], 0.5 on the diagonal / no data
% M.p              GLM estimate (their "x2_pMat_glm"): per target neuron i, the pair
%                  perturbations are shrunk towards the neuron's mean (precision-weighted
%                  posterior mean Q, eq. 3) and fit with a binomial probit GLM; then
%                  p(i,j) / (p(i,j) + p(j,i)) and clamped to [0.01 0.99]
% M.W              log M.p - log M.p' : log-likelihood ratio that i fires after j
%                  (antisymmetric), used by seq_likelihood
% Pairs are counted over all spikes in a window (not only neighbouring spikes), as in
% the original. Counting uses cumulative spike counts instead of their double loop.

n_after = zeros(n);
for w = 1:numel(windows)
    s = windows{w};
    if size(s, 1) < 2; continue; end
    [Hprev, X] = earlier_counts(s(:, 1), s(:, 2), n);
    n_after = n_after + X' * Hprev;   % (i, j): spikes of i preceded by spikes of j
end
n_pairs = n_after + n_after';
n_pairs(1:n + 1:end) = 0;

p_raw = n_after ./ n_pairs;
p_raw(isnan(p_raw)) = 0.5;
p_raw(1:n + 1:end) = 0.5;
p_raw = min(max(p_raw, 0.01), 0.99);

% GLM-based estimate (their code, row = target neuron i)
p_glm = 0.5 * ones(n);
for i = 1:n
    tot = sum(n_pairs(i, :));
    if tot == 0; continue; end
    m0 = (2 * sum(n_after(i, :)) - tot) / tot;
    sig0 = (sum(n_after(i, :)) * (1 - m0)^2 + (tot - sum(n_after(i, :))) * (-1 - m0)^2) / tot;
    j = find(n_pairs(i, :) > 0);
    q = (2 * n_after(i, j) - n_pairs(i, j)) ./ n_pairs(i, j);   % sample mean perturbation
    sigj = sum(n_pairs(i, j) .* (q - m0).^2) / tot;
    if numel(j) < 2 || sigj == 0 || sig0 == 0
        p_glm(i, j) = n_after(i, j) ./ n_pairs(i, j);
        continue
    end
    Q = (n_pairs(i, j) .* q / sigj + m0 / sig0) ./ (n_pairs(i, j) / sigj + 1 / sig0);
    yn = [n_after(i, j); n_pairs(i, j)]';
    warning('off', 'stats:glmfit:IterationLimit');
    warning('off', 'stats:glmfit:PerfectSeparation');
    b = glmfit(Q', yn, 'binomial', 'link', 'probit');
    p_glm(i, j) = glmval(b, Q', 'probit');
end
p = p_glm ./ (p_glm + p_glm');
p(isnan(p)) = 0.5;
p(1:n + 1:end) = 0.5;
p = min(max(p, 0.01), 0.99);

M.n_after = n_after;
M.n_pairs = n_pairs;
M.p_raw = p_raw;
M.p = p;
M.W = log(p) - log(p');
M.n_windows = numel(windows);
end


function [Hprev, X] = earlier_counts(t, id, n)
% Hprev(k, j): number of spikes of neuron j strictly earlier than spike k
% (spikes at the same time are not counted as earlier). X: one-hot (spikes x n).
N = numel(t);
[~, ~, g] = unique(t);
X = sparse(1:N, id, 1, N, n);
G = sparse(g, 1:N, 1, max(g), N) * X;                     % per distinct time
prev = [zeros(1, n); cumsum(full(G(1:end - 1, :)), 1)];
Hprev = prev(g, :);
end
