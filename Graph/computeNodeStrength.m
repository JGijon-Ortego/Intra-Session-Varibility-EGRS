function strength = computeNodeStrength(adjMatrix)
%COMPUTENODESTRENGTH Computes node strength (weighted degree) of each node.
% Node strength is the sum of the weights of the links connected to the
% node. Equivalent to strengths_und of the Brain Connectivity Toolbox.
%
% Inputs:
%   adjMatrix: NxN undirected weighted connection matrix (N = number of
%              nodes), with the weight of each connection m-n at position
%              (m,n). It should be symmetric. If the diagonal is not zero,
%              the diagonal values are included in the sum.
%
% Outputs:
%   strength: 1xN vector with the strength of each node.

strength = sum(adjMatrix);
end
