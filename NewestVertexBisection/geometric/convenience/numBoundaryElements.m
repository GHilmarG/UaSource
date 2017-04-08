function nBE = numBoundaryElements(mesh)
nBE = zeros(1,numBoundaries(mesh));
for j = 1:numBoundaries(mesh)
    nBE(j) = size(mesh.bd(j).elements,1);
end