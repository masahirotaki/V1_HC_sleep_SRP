function [pred, score] = apply_decoder(m, X)
% APPLY_DECODER  Predict classes of spike-count samples X with a fit_decoder model.
% score [nSamples x K]: Poisson log-likelihood ('bayes'), count-conditioned
% multinomial log-likelihood ('multinomial'), discriminant ('lda', 'pca_lda'),
% predicted class indicator ('pls') or SVM score ('svm').
% For 'bayes', X must use the same bin width as the training counts; 'multinomial'
% does not depend on bin width or overall gain.
% Samples whose scores tie across classes (e.g. no spikes, for 'multinomial') get a
% random class among the tied ones, so they count at chance.

if strcmp(m.method, 'bayes')
    score = X * log(m.lambda') - sum(m.lambda, 2)';
elseif strcmp(m.method, 'multinomial')
    score = X * log(m.P');
else
    Z = (X - m.mu) ./ m.sd;
    switch m.method
        case 'lda'
            score = Z * m.lda.W + m.lda.b;
        case 'pca_lda'
            score = (Z * m.coeff) * m.lda.W + m.lda.b;
        case 'pls'
            score = [ones(size(Z, 1), 1) Z] * m.beta;
        case 'svm'
            [~, score] = predict(m.svm, Z);
    end
end
[best, pred] = max(score, [], 2);
tied = find(sum(score >= best - 1e-12, 2) > 1);
for i = tied'
    cand = find(score(i, :) >= best(i) - 1e-12);
    pred(i) = cand(randi(numel(cand)));
end
end
