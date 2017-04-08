function [face2nodes,element2faces,varargout] = getHyperfaces(mesh)
%GETHYPERFACES    Returns all the hyperfaces of the mesh.
%   [FACE2NODES, ELEMENT2FACES, BOUNDARYONE2FACE, ...] = getHyperfaces(MESH)
%   yields the faces FACE2NODES and the following
%   indices with respect to FACE2NODES: 
%   ELEMENT2FACES(i,:) the indices of the hyperfaces of the i-th element,
%   whereas the hyperface ELEMENT2FACES(i,j) is the hyperface
%      ELEMENTS(i,(1:n)~=j).
%   This ordering is derived from the function SIMPLEXHYPERFACES.
%   BOUNDARYONE2FACE(i) the index of the i-th boundary element of bd(1)
%   and the same for the consecuting boundaries.
%
%   Works for arbitrary-dimensional meshes.
%
%   Uniform orientation of face2nodes is NOT guaranteed. The ordering
%   of the nodes in FACE2NODES is ascending so the orientations of the
%   faces of each element alternate.
%
%   If you want more flexibility use the function GETSUBSIMPLICES.
%
%   See also:
%   GETSUBSIMPLICES
%   GENMESH
%   GETBOUNDARY
%   GETEDGES
%
% Author: Josef Kemetmueller - 20.03.13

% This ordering is essential for the later equality to hold:
% element2faces(j,k) does not contain the node elements(j,k).
% so elements(j,k) \mcommentfont $\notin$ face2nodes(element2faces(j,k),:)
hyperfaces = simplexHyperfaces(dimMesh(mesh));
boundary = 1:dimMesh(mesh);
%% Compute the facedata
varargout = cell(1,numBoundaries(mesh));
[face2nodes, element2faces, varargout{:}] = ...
                    getSubsimplices(mesh, hyperfaces, boundary);


