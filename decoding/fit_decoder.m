function m = fit_decoder(X, y, K, method)
% FIT_DECODER  Train a population decoder on spike counts (apply with apply_decoder).
%
% X [nSamples x nNeurons] spike counts, y [nSamples x 1] labels 1..K
% method:
%   'bayes'    Poisson naive Bayes (memoryless Bayesian decoder, flat prior), raw counts.
%              m.lambda = expected count per class x neuron (rates x bin width).
%   'multinomial'  Bayesian decoder conditioned on the total spike count of the sample:
%              only WHICH cells fired is used, not how many spikes in total, so a global
%              gain change drops out. m.P(k, i) = share of spikes from cell i in class k,
%              P = (lambda + 1e-3) normalised over cells; log posterior = counts * log(P)'.
%   'lda'      linear discriminant, Ledoit-Wolf shrinkage covariance, z-scored counts
%   'pca_lda'  PCA (components explaining 90% variance, max 20) then shrinkage LDA
%   'pls'      PLS-DA: PLS regression onto one-hot classes (max 10 components), argmax
%   'svm'      linear SVM, one-vs-all (fitcecoc + templateLinear)
% z-scoring, PCA and PLS are fit on the training data only.

m.method = method;
m.K = K;
if any(strcmp(method, {'bayes', 'multinomial'}))
    lam = zeros(K, size(X, 2));
    for k = 1:K
        lam(k, :) = mean(X(y == k, :), 1);
    end
    if strcmp(method, 'bayes')
        m.lambda = max(lam, 1e-3);
    else
        m.P = (lam + 1e-3) ./ sum(lam + 1e-3, 2);
    end
    return
end

m.mu = mean(X, 1);
m.sd = std(X, 0, 1);
m.sd(m.sd == 0) = 1;
Z = (X - m.mu) ./ m.sd;
switch method
    case 'lda'
        m.lda = fit_shrinkage_lda(Z, y, K);
    case 'pca_lda'
        [coeff, ~, latent] = pca(Z);
        n_pc = min([find(cumsum(latent) / sum(latent) >= 0.9, 1), 20, size(coeff, 2)]);
        m.coeff = coeff(:, 1:n_pc);
        m.lda = fit_shrinkage_lda(Z * m.coeff, y, K);
    case 'pls'
        Y = double(y(:) == 1:K);
        n_comp = min([10, size(Z, 2), size(Z, 1) - 1]);
        [~, ~, ~, ~, m.beta] = plsregress(Z, Y, n_comp);
    case 'svm'
        m.svm = fitcecoc(Z, y, 'Learners', templateLinear('Learner', 'svm'), 'Coding', 'onevsall', ...
            'ClassNames', 1:K);
    otherwise
        error('Unknown decoder %s', method);
end
end


function lda = fit_shrinkage_lda(Z, y, K)
% Pooled within-class covariance shrunk towards a scaled identity (Ledoit & Wolf 2004).
p = size(Z, 2);
mu = zeros(K, p);
for k = 1:K
    mu(k, :) = mean(Z(y == k, :), 1);
end
Xc = Z - mu(y, :);
n = size(Xc, 1);
S = (Xc' * Xc) / n;
target = trace(S) / p;
d2 = norm(S - target * eye(p), 'fro')^2;
b2 = (sum(sum(Xc.^2, 2).^2) - n * norm(S, 'fro')^2) / n^2;
shrink = max(0, min(1, min(b2, d2) / max(d2, eps)));
Sigma = shrink * target * eye(p) + (1 - shrink) * S;
lda.W = Sigma \ mu';
lda.b = -0.5 * sum(mu' .* lda.W, 1);
lda.shrinkage = shrink;
end
