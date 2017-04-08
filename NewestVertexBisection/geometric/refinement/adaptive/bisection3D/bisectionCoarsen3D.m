function mesh = bisectionCoarsen3D(mesh, markedElements, mode)
%BISECTIONCOARSEN3D    Coarsens a bisection mesh
%   MESH = BISECTIONCOARSEN3D(MESH, markedElements) coarsens the elements
%   markedElements.
%   Default behaviour is to only coarsen edges, whose entire edge patches
%   have been marked for coarsening.
%   You can use 
%   MESH = BISECTIONCOARSEN3D(MESH, markedElements, 'force') to coarsen an
%   edge patch, even though only one of its elements has been marked.
%   If you want to mark all elements cou can use
%   MESH = BISECTIONCOARSEN3D(MESH, 'all');
%   Repeated application of above function will finally yield the original
%   mesh.
%
%   Works for 3-dimensional meshes in R^3 only.
%
%   Use the function BISECTIONCOARSEN instead of this function directly.
%
%   See also:
%   BISECTIONCOARSEN
%
% Author: Josef Kemetmueller - 20.03.13

%% Remark:
% The implementation relies on the fact that bisectionRefine3D
% always yields a mesh with the left child of an element being ordered in
% front of the right one and in this way keeping the implicit binary tree
% intact.
%%
nC = numCoordinates(mesh);
nE = numElements(mesh);
nB = numBoundaries(mesh);
if ischar(markedElements)
    assert(strcmpi(markedElements,'all'), 'Unknown option %s', markedElements);
    markedElements = 1:nE;
end
assert(isfield(mesh,'elementgeneration'), 'Mesh was not generated with genBisectionMesh. ');
%% Which nodes should be deleted?
% How many mesh.elements use the node?
node2numTet = accumarray(reshape(mesh.elements,[],1),1,[nC 1]);
if exist('mode','var')
    assert(strcmpi(mode,'force'), 'Unknown option %s', mode);
    isNodeOfMarkedEl = logical(accumarray(reshape(mesh.elements(markedElements,:),[],1),1,[nC 1]));
    newestnode2numTet = accumarray(mesh.elements(:,2),1,[nC 1]);
    % A vertex of a marked element can be removed, if it is the newest
    % vertex of all the mesh.elements of its nodepatch.
    nodesToRemove = isNodeOfMarkedEl & (node2numTet == newestnode2numTet);
else
    % How many marked mesh.elements use the node as their newest vertex?
    markednewestVert2numTet = accumarray(mesh.elements(markedElements,2),1,[nC 1]);
    % A vertex can be removed, if it is the newest vertex of all the
    % mesh.elements of its nodepatch.
    nodesToRemove = (node2numTet == markednewestVert2numTet);
end
%% Which elements are being coarsened?
element2neighbors = getNeighbors(mesh);
leftHalf = ((mesh.elementgeneration > 2) & nodesToRemove(mesh.elements(:,2)) ...
                                         & (element2neighbors(:,1) > (1:nE)'));
% Do not remove nodes whose elementgeneration is <=2.
nodesToRemove = false(nC,1);
nodesToRemove(mesh.elements(leftHalf,2)) = 1;
% Generate new numbers.
mesh.coordinates = mesh.coordinates(~nodesToRemove,:);
coordinates2newCoordinates = zeros(1,nC);
coordinates2newCoordinates(~nodesToRemove) = 1:numCoordinates(mesh);
%% Coarsen the triangulation.
rightHalf = element2neighbors(leftHalf,1); 
mesh.elements(leftHalf,:) = [mesh.elements(leftHalf,[1,3,4]),mesh.elements(rightHalf,1)];
mesh.elements(rightHalf,:) = [];
mesh.elementgeneration(leftHalf) = mesh.elementgeneration(leftHalf)-1;
mesh.elementgeneration(rightHalf) = [];
mesh.elements = coordinates2newCoordinates(mesh.elements);

%% Coarsen the boundary faces
%   There are only two possibilities:
%   Four triangles being melted into two trianges.
%   Two triangles being melted into one triangle
nBE = numBoundaryElements(mesh);
for j = 1:nB
    if ~isempty(mesh.bd(j).elements)        
        boundary2neighbors = getNeighbors(meshBd(mesh,j));
        leftHalf = (nodesToRemove(mesh.bd(j).elements(:,2)) & (boundary2neighbors(:,1) > (1:nBE(j))'));
        rightHalf = boundary2neighbors(leftHalf,1);
        
        mesh.bd(j).elements(leftHalf,:) = [mesh.bd(j).elements(leftHalf,[1,3]), ...
                                           mesh.bd(j).elements(rightHalf,1)];
        mesh.bd(j).elements(rightHalf,:) = [];
        mesh.bd(j).elements = coordinates2newCoordinates(mesh.bd(j).elements);
    end
end