function [edge2nodes, element2edges, varargout] = getEdges(mesh)
%GETEDGES    Returns all the edges of the mesh
%   [EDGE2NODES, ELEMENT2EDGES, BOUNDARYONE2EDGES, ...] =
%   GETEDGES(MESH) yields the edges EDGE2NODES and the following
%   indices with respect to EDGE2NODES:
%   ELEMENT2EDGES(i,:) the indices of the edges of the i-th element.
%   The element [1,2,3,4] yields the following edges in ascending
%   ordering lexicographically:
%   (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).
%   The ordering of these edges is derived from the function SIMPLEXEDGES.
%   The arrays BOUNDARYX2EDGES of a hyperface [1,2,3] is sorted in the same
%   lexicographical manner (1,2), (1,3), (2,3).
%   The ordering of the nodes of each edge of EDGE2NODES is ascending.
%
%   Works for arbitrary-dimensional meshes. 
%
%   If you want more flexibility use the function GETSUBSIMPLICES.
%
%   See also:
%   GETSUBSIMPLICES
%   GENMESH
%   GETBOUNDARY
%   GETHYPERFACES
%
% Author: Josef Kemetmueller - 20.03.13

%% Ordering of the edges
edgeorder = simplexEdges(dimMesh(mesh));
boundaryedgeorder = simplexEdges(dimMesh(mesh)-1);
%% Compute the edgedata.
varargout = cell(1,numBoundaries(mesh));
[edge2nodes, element2edges, varargout{:}] = ...
                       getSubsimplices(mesh, edgeorder, boundaryedgeorder);