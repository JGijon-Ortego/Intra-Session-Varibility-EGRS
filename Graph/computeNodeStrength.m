function strength = computeNodeStrength(adjMatrix)
%COMPUTENODESTRENGTH Computes node strength (weighted degree) of each node.
% Valid for undirected graphs without self-loops.
%
% Inputs:
%   adjMatrix: NxN matrix (N = number of nodes) containing the weight of
%              each connection m-n at position (m,n).
%
% Outputs:
%   strength: vector with the strength of each node.

sz = size(adjMatrix);

% Symmetrize and remove self-loops
idx = find(triu(squeeze(adjMatrix), 1));
mask = zeros(sz(1), sz(2));
mask(idx) = 1;
adjMatrix = adjMatrix .* mask;   % upper triangular matrix
adjMatrix = adjMatrix + adjMatrix';

% Strength of each node
strength = nansum(adjMatrix) / 2;
end
