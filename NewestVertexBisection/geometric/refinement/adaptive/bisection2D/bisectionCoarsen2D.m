function mesh = bisectionCoarsen2D(mesh, markedElements, mode)
%BISECTIONCOARSEN2D    Coarsens a bisection mesh
%   MESH = BISECTIONCOARSEN2D(MESH, markedElements) coarsens the elements
%   markedElements.
%   Default behaviour is to only coarsen edges, whose entire edge patches
%   have been marked for coarsening.
%   You can use 
%   MESH = BISECTIONCOARSEN2D(MESH, markedElements, 'force') to coarsen an
%   edge patch, even though only one of its elements has been marked.
%   If you want to mark all elements cou can use
%   MESH = BISECTIONCOARSEN2D(MESH, 'all');
%   Repeated application of above function will not always yield the original
%   mesh. Some initial triangulations that have the refinement edge chosen
%   in a special way will not work.
%   Here is an example of a mesh, that can't be fully coarsened after
%   refinement:
%       N = 5;
%       [X,Y] = pol2cart(linspace(0,2*pi*(1-1/N),N),1);
%       coordinates = [X(:), Y(:); 0,0];
%       elements = [(N+1)*ones(N,1), (1:N)', [(2:N)';1]];
%       mesh = genMesh(elements, coordinates, 'outwards');
%   If you prefer a mesh, that can be coarsened to its original state, use
%   the command MESH = GENBISECTIONMESH(MESH, 'forceRefine'); to generate
%   an initial mesh, that has this property.
%
%   Works for 2-dimensional meshes in R^n.
%
%   Use the function BISECTIONCOARSEN instead of this function directly.
%
%   See also:
%   BISECTIONCOARSEN
%
% Author: Josef Kemetmueller - 20.03.13

%% Remark:
% The implementation relies on the fact that bisectionRefine2D
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
assert(isfield(mesh,'N0'), 'Mesh was not generated with genBisectionMesh2D. ');
%% Which nodes should be deleted?
% How many mesh.elements use the node?
node2numEl = accumarray(reshape(mesh.elements,[],1),1,[nC 1]);
if exist('mode','var')
    assert(strcmpi(mode,'force'), 'Unknown option %s', mode);
    isNodeOfMarkedEl = logical(accumarray(reshape(mesh.elements(markedElements,:),[],1),1,[nC 1]));
    newestnode2numTet = accumarray(mesh.elements(:,2),1,[nC 1]);
    % A vertex of a marked element can be removed, if it is the newest
    % vertex of all the mesh.elements of its nodepatch.
    nodesToRemove = isNodeOfMarkedEl & (node2numEl == newestnode2numTet);
    nodesToRemove(1:mesh.N0) = false;
else
    % How many marked mesh.elements use the node as their newest vertex?
    markednewestVert2numEl = accumarray(mesh.elements(markedElements,2),1,[nC 1]);
    % A vertex can be removed, if it is the newest vertex of all the
    % mesh.elements of its nodepatch.
    nodesToRemove = (node2numEl == markednewestVert2numEl);
    nodesToRemove(1:mesh.N0) = false;
end
%% Which elements are being coarsened?
element2neighbors = getNeighbors(mesh);
leftHalf = (nodesToRemove(mesh.elements(:,2)) & (element2neighbors(:,1) > (1:nE)'));
nodesToRemove = false(nC,1);
nodesToRemove(mesh.elements(leftHalf,2)) = true;
% Generate new numbers.
mesh.coordinates = mesh.coordinates(~nodesToRemove,:);
coordinates2newCoordinates = zeros(1,nC);
coordinates2newCoordinates(~nodesToRemove) = 1:numCoordinates(mesh);
%% Coarsen the triangulation.
rightHalf = element2neighbors(leftHalf,1); 
mesh.elements(leftHalf,:) = [mesh.elements(leftHalf,[1,3]),mesh.elements(rightHalf,1)];
mesh.elements(rightHalf,:) = [];
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
        
        mesh.bd(j).elements(leftHalf,:) = [mesh.bd(j).elements(leftHalf, 1), ...
                                           mesh.bd(j).elements(rightHalf,2)];
        mesh.bd(j).elements(rightHalf,:) = [];
        mesh.bd(j).elements = coordinates2newCoordinates(mesh.bd(j).elements);
    end
end