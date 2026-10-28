function [clcAverage, clcNodal] = computeClusteringCoefficient(W)
%COMPUTECLUSTERINGCOEFFICIENT Computes the clustering coefficient of the
% graph in W (valid for directed graphs).
% Portions of this function are derived from the Brain Connectivity
% Toolbox (BCT), Copyright (c) 2019, BCT Contributors, distributed under
% the BSD 3-Clause License (see THIRD_PARTY_LICENSES.md).
% Reference: Rubinov M, Sporns O (2010) NeuroImage 52:1059-1069.
%
% The following definition of the clustering coefficient is used:
%
%        ___   ___
%        \     \    Wik * Wli * Wkl
%        /__   /__
%        k~=i (l~=k & l~=i)
%   Ci = -----------------------------
%            ___   ___
%            \     \    Wik * Wil
%            /__   /__
%            k~=i (l~=k & l~=i)
%
% Stam CJ et al. (2009). Graph theoretical analysis of magnetoencephalographic
% functional connectivity in Alzheimer's disease. Brain, 132:213-224.
%
% Inputs:
%   W: NxN matrix (N = number of nodes) with the weight of each
%      connection m-n at position (m,n).
%
% Outputs:
%   clcAverage: global (mean) clustering coefficient.
%   clcNodal:   clustering coefficient of each node.
%
% For the global measure, surrogate testing should be applied afterwards.

A = W ~= 0;                               % adjacency matrix
S = W.^(1/3) + (W.').^(1/3);              % symmetrized weights matrix ^1/3
K = sum(A + A.', 2);                      % total degree (in + out)
cyc3 = diag(S^3) / 2;                     % number of 3-cycles (i.e., directed triangles)
K(cyc3 == 0) = inf;                       % if no 3-cycles exist, make C = 0 (via K = inf)
CYC3 = K .* (K - 1) - 2 * diag(A^2);      % number of all possible 3-cycles
clcNodal = cyc3 ./ CYC3;                  % clustering coefficient (nodal)
clcAverage = mean(clcNodal);              % clustering coefficient (global)
end
