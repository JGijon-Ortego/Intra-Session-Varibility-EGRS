function surrogate = surrogateGraph(adjMatrix)
%SURROGATEGRAPH Generates a surrogate graph by shuffling the connections.
% Only a single shuffle is performed; it does not preserve the node
% degree distribution.
%
% Inputs:
%   adjMatrix: NxN matrix (N = number of nodes) with the weight of each
%              connection m-n at position (m,n). The graph must be
%              undirected (symmetric matrix).
%
% Outputs:
%   surrogate: surrogate graph (symmetric NxN matrix).

surrogateUpper = zeros(size(adjMatrix));
idxUpper = find(triu(adjMatrix, 1));
valuesUpper = adjMatrix(idxUpper);
valuesShuffled = valuesUpper(randperm(length(valuesUpper)));
surrogateUpper(idxUpper) = valuesShuffled;
surrogate = surrogateUpper + surrogateUpper';
end
