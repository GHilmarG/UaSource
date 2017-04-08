function mesh = redRefine3D(mesh)
%REDREFINE3D    Refines a given mesh uniformly using red/regular refinement
%   MESH = REDREFINE3D(MESH) executes one uniform refinement step.
%
%   The process results in the following eight tetrahedra per element
%   (v1,  v12, v13, v14),     (v12, v13, v14, v24),
%   (v12, v2,  v23, v24),     (v12, v13, v23, v24),
%   (v13, v23, v3,  v34),     (v13, v23, v24, v34),
%   (v14, v24, v34, v4 ) and  (v13, v14, v24, v34).
%   vi being the i-th vertex of the tetrahedron and
%   vij being the midpoint of the edge (i,j).
%
%   Works for 3-dimensional meshes in R^3 only.
%
% Author: Josef Kemetmueller - 20.03.13

nB = numBoundaries(mesh);
nC = numCoordinates(mesh);
[edge2nodes,element2edges,boundary2edges{1:nB}] = getEdges(mesh);
% Generate edge midpoints
midpoints = @(C,edges) 0.5*(C(edges(:,1),:)+C(edges(:,2),:));
mesh.coordinates = [mesh.coordinates; midpoints(mesh.coordinates,edge2nodes)];
% Sorting from getEdges: [1,2; 1,3; 1,4; 2,3; 2,4; 3,4]
newElements = @(el, el2edgecenter) ...
                [el(:,1), el2edgecenter(:,[1,2,3]);...
                 el2edgecenter(:, 1),el(:,2),el2edgecenter(:,[4,5]);...
                 el2edgecenter(:,[2,4]),el(:,3),el2edgecenter(:, 6 );...
                 el2edgecenter(:,[3,5,6]),el(:,4);...
                 el2edgecenter(:,[1,2,3,5]);...
                 el2edgecenter(:,[1,2,4,5]);...
                 el2edgecenter(:,[2,4,5,6]);...
                 el2edgecenter(:,[2,3,5,6]);...
                ];
mesh.elements = newElements(mesh.elements, nC+element2edges);
        
% Refine boundary faces
for j=1:nB
    bd = mesh.bd(j).elements;
    if ~isempty(bd)
        mesh.bd(j).elements = [bd(:,1),nC+boundary2edges{j}(:,[1,2]);...
                               bd(:,2),nC+boundary2edges{j}(:,[3,1]);...
                               bd(:,3),nC+boundary2edges{j}(:,[2,3]);...
                               nC+boundary2edges{j}(:,[3,2,1])];
    else
        mesh.bd(j).elements = [];
    end
end
% Remove elementgeneration field to mark the mesh as not being bisectable
% anymore
if isfield(mesh, 'elementgeneration')
    mesh = rmfield(mesh,'elementgeneration');
end