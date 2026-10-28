function [lambda, efficiency, ecc, radius, diameter, D] = computePathLength(adjMatrix)
%COMPUTEPATHLENGTH Computes the path length and other integration measures
% of a graph. Requires the Brain Connectivity Toolbox (BCT) on the path.
% Portions of this function are derived from the Brain Connectivity
% Toolbox (BCT), Copyright (c) 2019, BCT Contributors, distributed under
% the BSD 3-Clause License (see THIRD_PARTY_LICENSES.md).
% Reference: Rubinov M, Sporns O (2010) NeuroImage 52:1059-1069.
%
% Inputs:
%   adjMatrix: NxN matrix (N = number of nodes) with the weight of each
%              connection m-n at position (m,n).
%
% Outputs:
%   lambda:     characteristic path length.
%   efficiency: global efficiency of the network.
%   ecc:        eccentricity of each node.
%   radius:     radius of the network.
%   diameter:   diameter of the network.
%   D:          matrix of shortest path distances.

L = 1 ./ adjMatrix;
L(find(1 ./ adjMatrix == Inf)) = 0;
[D, ~] = distance_wei(L);
[lambda, efficiency, ecc, radius, diameter] = charpath(D);
end
