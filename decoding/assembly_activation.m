function R = assembly_activation(Z, W)
% ASSEMBLY_ACTIVATION  Reactivation (activation) strength of each assembly per bin:
% R_k(t) = z(t)' P_k z(t) with P_k = w_k w_k' and its diagonal set to 0, so a neuron's own
% firing does not count, only co-firing of neurons with weights in the pattern
% (Lopes-dos-Santos et al. 2013).
%
% Z [nBins x nNeurons] z-scored counts; W [nNeurons x nAssemblies] (assembly_templates)
% R [nBins x nAssemblies]

R = (Z * W).^2 - (Z.^2) * (W.^2);
end
