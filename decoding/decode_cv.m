function pred = decode_cv(X, y, fold, method)
% DECODE_CV  Cross-validated class prediction from population spike counts.
%
%   pred = decode_cv(X, y, fold, method)
%
% X     [nSamples x nNeurons] spike counts
% y     [nSamples x 1] class labels 1..K
% fold  [nSamples x 1] test fold of each sample (make_folds); every fold is
%       predicted by a model trained on all other folds
% method  see fit_decoder ('bayes', 'lda', 'pca_lda', 'pls', 'svm')

K = max(y);
pred = zeros(size(y));
for f = unique(fold(:))'
    test = fold == f;
    m = fit_decoder(X(~test, :), y(~test), K, method);
    pred(test) = apply_decoder(m, X(test, :));
end
end
