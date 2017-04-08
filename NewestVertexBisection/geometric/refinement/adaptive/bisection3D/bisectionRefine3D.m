function mesh = bisectionRefine3D(mesh, markedElements)
%BISECTIONREFINE3D    Refines a mesh via bisection
%   MESH = BISECTIONREFINE3D(MESH, MARKEDELEMENTS) refines the elements
%   MARKEDELEMENTS via recursive bisection.
%
%   Not all meshes are compatibly divisible using this refinement scheme,
%   as it is required for an initial triangulation that for every two
%   neighboring elements T, T' there holds:
%      - If  one of the refinement edges of T or  T' lies in $T \cap T'$,
%        then T and T' must be reflected neighbors
%      - If none of the refinement edges of T and T' lies in $T \cap T'$,
%        their neighboring children must be reflected neighbors.
%   To obtain such an initial mesh from an arbitrary one, use the function
%   GENBISECTIONMESH3D.
%
%   The orientations of the refined elements and boundary faces will NOT be
%   uniform. This is a limitation of the underlying scheme. The
%   ordering/orientation of the elements and the nodes MUST NOT be changed,
%   as this will mess up the refinement routine and will probably lead to
%   MATLAB crashing. If you need a positively oriented mesh, you can get an
%   oriented copy of the mesh using ORIENTMESH.
%
%   Works for 3-dimensional meshes in R^3 only.
%
%   Use the function BISECTIONREFINE instead of this function directly.
%
%   See also:
%   BISECTIONREFINE
%   ORIENTMESH
%
% Author: Josef Kemetmueller - 20.03.13

nB = numBoundaries(mesh);
if islogical(markedElements)
    markedElements = find(markedElements);
end
if ischar(markedElements)
    assert(strcmpi(markedElements,'all'), 'Unknown option %s', markedElements);
    markedElements = 1:numElements(mesh);
end
assert(dimMesh(mesh)==3 && size(mesh.elementgeneration,2)==1 && dimSpace(mesh)==3, ...
       'Mesh dimensions not suitable for this algorithm.');
assert(~isempty(markedElements) && all(markedElements>=0) && all(markedElements<=numElements(mesh)),...
       'markedElements must be in range of number of elements.');
% Generate helping structures
[face2nodes,element2faces,boundary2face{1:nB}] = getHyperfaces(mesh);
element2neigh = getNeighbors(mesh);

assert(size(unique(element2faces(:)),1)==max(element2faces(:)), ...
       'The boundary faces don''t match the elements. Was the mesh not generated using genBisectionMesh?');
nF = size(face2nodes,1);
element2bdNumber = cell(1,nB);
nBE = numBoundaryElements(mesh);
for j = 1:nB
    bdNumber = zeros(nF,1);
    bdNumber(boundary2face{j}) = 1:nBE(j);
    element2bdNumber{j} = reshape(bdNumber(element2faces),[],4);
end
%*** Rufe MEX-Funktion
if (nB==0)
    [mesh.elements,mesh.elementgeneration,mesh.coordinates] = bisectionRefine3DRecursiveC(mesh.elements,mesh.elementgeneration,mesh.coordinates,element2neigh,markedElements);
else
    boundaries = struct2cell(mesh.bd);
    [mesh.elements,mesh.elementgeneration,mesh.coordinates,bd{1:nB}] = bisectionRefine3DRecursiveC(mesh.elements,mesh.elementgeneration,mesh.coordinates,element2neigh,boundaries{1:nB},element2bdNumber{1:nB},markedElements);
    mesh.bd = cell2struct(bd,'elements');
end