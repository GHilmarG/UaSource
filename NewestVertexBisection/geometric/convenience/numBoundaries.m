function nB = numBoundaries(mesh)
nB = 0;
if isfield(mesh,'bd')
    nB = length(mesh.bd);
end