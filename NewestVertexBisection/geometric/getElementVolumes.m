function volumes = getElementVolumes(mesh)
%GETELEMENTVOLUMES    Volumes(/areas/lengths) of the simplices in the mesh.
%   VOLUMES = GETELEMENTVOLUMES(MESH) returns a vector containing the
%   volumes(/areas/lengths) of each simplex.
%
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 16.12.2013

nE = numElements(mesh);

if(dimMesh(mesh)>dimSpace(mesh))
    warning(['Your %dD-simplices are lying in %dD-space.\n',...
             'The only correct ''volume'' is zero.'], dimMesh(mesh), dimSpace(mesh));
    volumes = zeros(nE,1);
    return;
end
if (dimMesh(mesh)==dimSpace(mesh))
    volumes = abs(getSignedElementVolumes(mesh));
else % Our mesh is a manifold in some higher dimensional space.
    X = @(d) mesh.coordinates(mesh.elements(:,d),:);
    lengths = @(X) sqrt(sum(X.^2,2));
    %% 
    if (dimMesh(mesh)==0) % counting measure
        volumes = ones(nE,1);
    elseif (dimMesh(mesh)==1) % Lengths are easy
        volumes = lengths(X(2)-X(1));
    elseif(dimMesh(mesh)==2 && dimSpace(mesh)==3) % Areas in 3D using cross product.
        volumes = lengths(cross(X(2)-X(1),X(3)-X(1),2))/2;
    else % Using general rule.
        volumes = zeros(nE,1);
        for i = 1:nE
            nodes = mesh.elements(i,:);
            MT = [ones(1,dimMesh(mesh)+1); mesh.coordinates(nodes,:)'];
            volumes(i) = 1/factorial(dimMesh(mesh))*sqrt(det(MT'*MT));
        end
    end
end