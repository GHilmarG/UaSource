function [subs2nodes, element2subs, varargout] = getSubsimplices(mesh, subsimplex, bdSubsimplex)
%GETSUBSIMPLICES    Returns subsimplices of a given pattern.
%   [SUBS2NODES, ELEMENT2SUBS, BOUNDARYONE2SUBS, ...] =
%   GETSUBSIMPLICES(MESH, SUBSIMPLEX, BDSUBSIMPLEX) yields the subsimplices
%   SUBS2NODES and the following indices with respect to SUBS2NODES:
%   ELEMENT2SUBS(i,:) the indices of the subsimplices of the i-th element,
%   in order of the rows of SUBSIMPLEX and BDSUBSIMPLEX respectively!
%   This means the call
%   [SUBS2NODES, ELEMENT2SUBS] = GETSUBSIMPLICES(mesh, [2,3; 1,3; 1,2])
%   will result in ELEMENT2SUBS(i,1) corresponding to the edge
%   ELEMENTS(i,[2,3]). Keep in mind though, that the ordering of the nodes
%   of each subsimplex of SUBS2NODES will ALWAYS be ascending. So the
%   orientation of SUBS2NODES will not be the same as in SUBSIMPLEX, as the
%   calls
%   subs2nodes = getSubsimplices(mesh, [1,2]) and
%   subs2nodes = getSubsimplices(mesh, [2,1]) will have the same results.
%
%   Works for arbitrary-dimensional meshes.
%
%   See also:
%   GETHYPERFACES
%   GETEDGES
%
% Author: Josef Kemetmueller - 20.03.13
nE = numElements(mesh);
nB = numBoundaries(mesh);
nBE = numBoundaryElements(mesh);
numSubs = size(subsimplex,1);
dimSubs = size(subsimplex,2);
assert( (max(subsimplex(:)) <= dimMesh(mesh)+1) & ...
        (min(subsimplex(:)) > 0), 'Subsimplex must be in range 1:dimMesh+1');

subsimplices = reshape(mesh.elements(:,subsimplex), nE*numSubs, dimSubs);

if exist('bdSubsimplex','var')
    numBdSubs = size(bdSubsimplex,1);
    dimBdSubs = size(bdSubsimplex,2);
    assert(dimSubs == dimBdSubs, 'bdSubsimplex must be of same dimension as subsimplex.');
    assert( (max(bdSubsimplex(:)) <= dimMesh(mesh)) & ...
            (min(subsimplex(:)) > 0), 'Subsimplex must be in range 1:dimMesh');
    bSubInd = cumsum([nE*numSubs, numBdSubs*nBE]);
    % Append the boundary faces
    for j = 1:nB
        boundary = mesh.bd(j).elements;
        if ~isempty(boundary)
            subsimplices(bSubInd(j)+1:bSubInd(j+1),:) = reshape(boundary(:,bdSubsimplex), [], dimSubs);
        end
    end
else
    numBdSubs = 0;
end
%% Do the actual work
[subs2nodes, ~, J] = unique(sort(subsimplices,2), 'rows');
element2subs = reshape(J(1:nE*numSubs),nE,numSubs);
%% Generate BOUNDARYX2SUBS
if numBdSubs > 0
    varargout = cell(1,nB);
    for j = 1:nB
        varargout{j} = reshape(J(bSubInd(j)+1:bSubInd(j+1)),nBE(j),numBdSubs);%boundary2subs{j}
    end
end