function A = assembly_templates(Z, method, seed)
% ASSEMBLY_TEMPLATES  Cell-assembly patterns from z-scored binned spike counts
% (Lopes-dos-Santos, Ribeiro & Tort 2013; Peyrache et al. 2009).
%
%   A = assembly_templates(Z)            % 'ica' (default)
%   A = assembly_templates(Z, 'pca')
%
% Z       [nBins x nNeurons] z-scored spike counts (each neuron mean 0, SD 1)
% Number of assemblies = eigenvalues of the correlation matrix above the Marchenko-Pastur
% bound lambda_max = (1 + sqrt(nNeurons / nBins))^2.
%   'pca'  patterns = the significant eigenvectors
%   'ica'  FastICA (symmetric, tanh) inside the significant PC subspace; patterns = columns
%          of the unmixing projected back to neuron space
% Patterns are unit length, signed so that the largest absolute weight is positive.
% A.W [nNeurons x nAssemblies], A.lambda (all eigenvalues), A.lambda_max, A.method
% Activation strength: assembly_activation(Z, A.W).

if nargin < 2; method = 'ica'; end
if nargin < 3; seed = 1; end
[nb, nn] = size(Z);
C = (Z' * Z) / nb;
[E, L] = eig((C + C') / 2);
[lam, o] = sort(diag(L), 'descend');
E = E(:, o);
lmax = (1 + sqrt(nn / nb))^2;
k = sum(lam > lmax);
A.lambda = lam;
A.lambda_max = lmax;
A.method = method;
if k == 0
    A.W = zeros(nn, 0);
    return
end
Ek = E(:, 1:k);
switch method
    case 'pca'
        W = Ek;
    case 'ica'
        Y = (Z * Ek) ./ sqrt(lam(1:k))';          % whitened projections (nBins x k)
        U = fastica_sym(Y, seed);                 % k x k, rows = unmixing vectors
        W = Ek * (U' ./ sqrt(lam(1:k)));          % back to neuron space
    otherwise
        error('Unknown method %s', method);
end
W = W ./ vecnorm(W);
[~, im] = max(abs(W), [], 1);
sg = sign(W(sub2ind(size(W), im, 1:size(W, 2))));
A.W = W .* sg;
end


function U = fastica_sym(Y, seed)
% symmetric FastICA with tanh contrast on whitened data Y (nBins x k)
rng(seed);
k = size(Y, 2);
U = orth(randn(k));
for it = 1:500
    S = Y * U';                                   % nBins x k
    g = tanh(S); gp = 1 - g.^2;
    Un = (g' * Y) / size(Y, 1) - mean(gp, 1)' .* U;
    [uu, ~, vv] = svd(Un, 'econ');
    Un = uu * vv';                                % symmetric decorrelation
    if max(abs(abs(diag(Un * U')) - 1)) < 1e-7; U = Un; break; end
    U = Un;
end
end
