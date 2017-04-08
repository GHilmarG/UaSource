function volumes = getSignedElementVolumes(mesh)
%GETSIGNEDELEMENTVOLUMES    Signed volumes(/areas/lengths) of the simplices
% in the given mesh.
%   VOLUMES = GETSIGNEDELEMENTVOLUMES(MESH) returns a vector containing the
%   signed volumes(/areas/lengths) of each simplex.
%
%   Works for n-dimensional meshes in R^n for arbitrary n.
%
%   Author: Josef Kemetmueller - 16.12.2013

nE = numElements(mesh);
assert(dimMesh(mesh)==dimSpace(mesh),['Your %dD-simplices are lying in %dD-space.',...
                                  'There is no signed volume.'],dimMesh(mesh),dimSpace(mesh));

X = cell(1,dimMesh(mesh)+1);
for d = 1:dimMesh(mesh)+1
    X{d} = mesh.coordinates(mesh.elements(:,d),:);
end
switch dimMesh(mesh)
    case 1 % Edge lengths are easy
        volumes = X{2}-X{1};
    case 2 % Manually compute 2D-determinants 
        d21 = X{2}-X{1};
        d31 = X{3}-X{1};
        volumes = 1/2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
    case 3 % The cross product can be used to compute the oriented volume.
        N{1} = -cross(X{3}-X{2},X{4}-X{2},2);
        volumes = (1/6)*dot(N{1},X{1}-X{2},2);
    otherwise % Generalized simplex volume.
        volumes = zeros(nE,1);
        for i = 1:nE
            nodes = mesh.elements(i,:);
            MT = [ones(1,dimMesh(mesh)+1); mesh.coordinates(nodes,:)'];
            volumes(i) = 1/factorial(dimMesh(mesh))*det(MT);
        end
end