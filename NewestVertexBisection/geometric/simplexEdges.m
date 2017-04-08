function edges = simplexEdges(dim)
%SIMPLEXEDGES    returns the edges of a simplex.
%   EDGES = SIMPLEXEDGES(DIM) returns the nchoosek(dim+1,2) edges of a
%   dim-dimensional simplex [1:dim+1]. The ordering of the nodes in an edge
%   is ascending and the edges are sorted lexicographically.
%
%   This function is meant to fix an ordering for the nodes of a simplex,
%   so that all functions depending on this ordering remain consistent. We
%   rely on this ordering for redRefine3D.
%
% Author: Josef Kemetmueller - 20.03.14

[I,J] = find(tril(ones(dim+1),-1));
edges = [J,I];