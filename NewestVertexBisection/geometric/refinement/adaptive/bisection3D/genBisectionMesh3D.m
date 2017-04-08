function Bmesh = genBisectionMesh3D(mesh)
%GENBISECTIONMESH3D    Generates a mesh, which is valid for adaptive bisection.
%   MESH = GENBISECTIONMESH3D(MESH) Returns a refined mesh which is valid for
%   adaptive bisection. 
%   For every tetrahedron it generates the following
%   twelve tetrahedra, using S as the element barycenter, and S_i as the
%   barycenter of the i-th face.
%
%   T_11 = (v_2, S_1, S, v_3)      T_21 = (v_1, S_2, S, v_3)
%   T_12 = (v_2, S_1, S, v_4)      T_22 = (v_1, S_2, S, v_4)
%   T_13 = (v_3, S_1, S, v_4)      T_23 = (v_3, S_2, S, v_4)
%
%   T_31 = (v_1, S_3, S, v_2)      T_41 = (v_1, S_4, S, v_2)
%   T_32 = (v_1, S_3, S, v_4)      T_42 = (v_1, S_4, S, v_3)
%   T_33 = (v_2, S_3, S, v_4)      T_43 = (v_2, S_4, S, v_3)
%
%   This yields a mesh that guarantees:
%     - for every two neighboring elements T, T' there holds:
%        - If  one of the refinement edges of T or  T' lies in $T \cap T'$,
%          then T and T' must be reflected neighbors.
%        - If none of the refinement edges of T and T' lies in $T \cap T'$,
%          their neighboring children must be reflected neighbors.
%   For more information on underlying mathematics see:
%       Rob Stevenson: The completion of locally refined simplicial
%       partitions created by bisection. Math. Comput. 77(261): 227-241
%       (2008)
%
%   Works for 3-dimensional meshes in R^3 only.
%
%   Use the function GENBISECTIONMESH instead of this function directly.
%
%   See also:
%   GENBISECTIONMESH
%
% Author: Josef Kemetmueller - 20.03.13

nC = numCoordinates(mesh);
nE = numElements(mesh);
nB = numBoundaries(mesh);
assert(dimSpace(mesh)==3, 'mesh must be embedded in 3D.');
assert(dimMesh(mesh)==3, 'mesh must be tetrahedral mesh.');
assert(~isfield(mesh,'elementgeneration'), ...
            ['Mesh seems to already be compatible. ',...
             '(Field ''elementgeneration'' already exists.)']);
% Get face information with built in correct ordering.
[face2nodes,element2faces,boundary2faces{1:nB}] = getHyperfaces(mesh);
% Ascending ordering of the nodes provided getHyperfaces is important!
nF = size(face2nodes,1);
% Generation of the barycenters and new coordinates.
elementcenters = getBarycenters(mesh);
facecenters = getBarycenters(genMesh(face2nodes,mesh.coordinates));
Bmesh.coordinates = [mesh.coordinates;elementcenters;facecenters];
% New node numbers
element2newNodes = nC+(1:nE)';
face2newNodes = nE+nC+(1:nF)';
% Generate the new elements according to the rule described above.
Bmesh.elements = [face2nodes(element2faces,1),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,2);...
                 face2nodes(element2faces,1),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,3);...
                 face2nodes(element2faces,2),face2newNodes(element2faces,:),repmat(element2newNodes,4,1),face2nodes(element2faces,3)];
Bmesh.elementgeneration = 2+zeros(12*nE,1);
% Generate the boundary faces so they match the simplices
for j = 1:nB
    if ~isempty(mesh.bd(j).elements)
        Bmesh.bd(j).elements = [face2nodes(boundary2faces{j},1),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},2);...
                               face2nodes(boundary2faces{j},1),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},3);...
                               face2nodes(boundary2faces{j},2),face2newNodes(boundary2faces{j}),face2nodes(boundary2faces{j},3)];
    else
        Bmesh.bd(j).elements = zeros(0,dimMesh(mesh));
    end
end