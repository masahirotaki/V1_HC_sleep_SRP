function R = decode_with_null(X, y, group, methods, opts)
% DECODE_WITH_NULL  K-fold cross-validated decoding plus a label-permutation null.
%
%   R = decode_with_null(X, y, group, methods, opts)
%
% X [nSamples x nNeurons] counts, y [nSamples x 1] labels 1..K, group [nSamples x 1]
% trial (sequence) of each sample. All samples of a trial are in the same test fold,
% so a decoder is never tested on a trial it saw in training.
% Null: labels are shuffled WITHIN each trial (every trial keeps one sample per
% label), and the whole cross-validation is repeated; p = (1 + #null >= acc) / (1 + n).
%
% opts.n_folds (10), .n_perm (200, scalar or one per method), .bin_width (s; for the
% decoded-time error of ordered labels, [] = not computed), .seed (1)
%
% R(m): .method .acc .chance .null_acc .p .pred .confusion (row = true, normalised)
%       .mae / .null_mae (mean |decoded - true| label distance x bin_width)

if ~isfield(opts, 'n_folds'); opts.n_folds = 10; end
if ~isfield(opts, 'n_perm'); opts.n_perm = 200; end
if ~isfield(opts, 'bin_width'); opts.bin_width = []; end
if ~isfield(opts, 'seed'); opts.seed = 1; end
n_perm = opts.n_perm .* ones(1, numel(methods));

rng(opts.seed);
[~, ~, g] = unique(group);
fold = make_folds(g, opts.n_folds);
members = accumarray(g, (1:numel(g))', [], @(i) {i});
perms = cell(max(n_perm), 1);
for i = 1:max(n_perm)
    idx = (1:numel(y))';
    for k = 1:numel(members)
        idx(members{k}) = members{k}(randperm(numel(members{k})));
    end
    perms{i} = idx;
end

K = max(y);
for m = 1:numel(methods)
    t0 = tic;
    pred = decode_cv(X, y, fold, methods{m});
    null_acc = nan(n_perm(m), 1);
    null_mae = nan(n_perm(m), 1);
    meth = methods{m};
    parfor i = 1:n_perm(m)
        yp = y(perms{i}); %#ok<PFBNS>
        pp = decode_cv(X, yp, fold, meth);
        null_acc(i) = mean(pp == yp);
        null_mae(i) = mean(abs(pp - yp));
    end
    R(m).method = meth; %#ok<AGROW>
    R(m).acc = mean(pred == y);
    R(m).chance = 1 / K;
    R(m).null_acc = null_acc;
    R(m).p = (1 + sum(null_acc >= R(m).acc)) / (1 + n_perm(m));
    R(m).pred = pred;
    R(m).confusion = confusionmat(y, pred, 'Order', 1:K) ./ accumarray(y, 1, [K 1]);
    if ~isempty(opts.bin_width)
        R(m).mae = mean(abs(pred - y)) * opts.bin_width;
        R(m).null_mae = null_mae * opts.bin_width;
        R(m).p_mae = (1 + sum(R(m).null_mae <= R(m).mae)) / (1 + n_perm(m));
    end
    fprintf('   %-8s acc %.3f (chance %.3f, null 95%% %.3f, p = %.3g)  [%.0f s]\n', meth, R(m).acc, ...
        1 / K, prctile(null_acc, 95), R(m).p, toc(t0));
end
end


function fold = make_folds(g, n_folds)
% Random, balanced assignment of whole trials to folds.
n_groups = max(g);
gf = mod(randperm(n_groups) - 1, n_folds) + 1;
fold = gf(g)';
end
