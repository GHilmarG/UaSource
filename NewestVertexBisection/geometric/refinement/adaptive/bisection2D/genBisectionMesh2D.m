function Bmesh = genBisectionMesh2D(mesh)
%GENBISECTIONMESH2D    Generates a mesh, which is valid for adaptive bisection.
%   MESH = GENBISECTIONMESH2D(MESH) Returns a refined mesh which is valid for
%   adaptive bisection. This is only necessary if you want the mesh to be
%   able to be coarsened to a full extent. 
%
%   For every triangle (v_1, v_2, v_3) it generates the following three
%   triangles, using S as the element barycenter
%           T_1 = (v_1, S, v_2)
%           T_2 = (v_2, S, v_3)
%           T_3 = (v_3, S, v_1)
%
%   Works for 2-dimensional meshes in R^n.
%
%   Use the function GENBISECTIONMESH instead of this function directly.
%
%   See also:
%   GENBISECTIONMESH
%
% Author: Josef Kemetmueller - 20.03.13
nC = numCoordinates(mesh);
nE = numElements(mesh);
assert(dimMesh(mesh)==2, 'mesh must be 2D.');
% New coordinates
Bmesh.coordinates = [mesh.coordinates; getBarycenters(mesh)];
Bmesh.N0 = size(Bmesh.coordinates,1);
% New node numbers
element2newNodes = nC+(1:nE)';
% Generate the new elements according to the rule described above.
Bmesh.elements = [mesh.elements(:,1),element2newNodes,mesh.elements(:,2); ...
                  mesh.elements(:,2),element2newNodes,mesh.elements(:,3); ...
                  mesh.elements(:,3),element2newNodes,mesh.elements(:,1)];
% No refinement needed for boundaries
Bmesh.bd = mesh.bd;