function mesh = bisectionCoarsen1D(mesh, markedElements)
%BISECTIONCOARSEN1D    Coarsens a bisection mesh
%   MESH = BISECTIONCOARSEN1D(MESH, markedElements) coarsens the elements
%   markedElements if its neighbor is also marked and of the same
%   generation.
%   MESH = BISECTIONCOARSEN1D(MESH, 'all');
%
%   Works for 1-dimensional meshes in R^n.
%
%   Use the function BISECTIONCOARSEN instead of this function directly.
%
%   See also:
%   BISECTIONCOARSEN
%
% Author: Josef Kemetmueller - 20.03.13

%% Remark:
% The implementation relies on the fact that bisectionRefine1D
% always yields a mesh with the left child of an element being ordered in
% front of the right one and in this way keeping the implicit binary tree
% intact. 
%%
nC = numCoordinates(mesh);
nE = numElements(mesh);
if ischar(markedElements)
    assert(strcmpi(markedElements,'all'), 'Unknown option %s', markedElements);
    markedElements = 1:nE;
end
assert(isfield(mesh,'N0'), 'Mesh was not generated with genBisectionMesh1D. ');
%% Which nodes should be deleted?
newest = @(E)  accumarray(E(:,2), E(:,2)>E(:,1), [nC 1]) & ...
               accumarray(E(:,1), E(:,1)>E(:,2), [nC 1]);
nodesToRemove = newest(mesh.elements(markedElements,:));
nodesToRemove(1:mesh.N0) = false;
%% Which elements are being coarsened?
element2neighbors = getNeighbors(mesh);
leftHalf = (nodesToRemove(mesh.elements(:,2)) & (element2neighbors(:,1) > (1:nE)'));
nodesToRemove = false(nC,1);
nodesToRemove(mesh.elements(leftHalf,2)) = true;
%% Generate new numbers.
mesh.coordinates = mesh.coordinates(~nodesToRemove,:);
coordinates2newCoordinates = zeros(1,nC);
coordinates2newCoordinates(~nodesToRemove) = 1:numCoordinates(mesh);
%% Coarsen the edges.
rightHalf = element2neighbors(leftHalf,1); 
mesh.elements(leftHalf,:) = [mesh.elements(leftHalf,1),mesh.elements(rightHalf,2)];
mesh.elements(rightHalf,:) = [];
mesh.elements = coordinates2newCoordinates(mesh.elements);
%% 1D boundaries don't change! YIPPIEH!