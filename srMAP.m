function Znew = srMAP(xNew,yNew,F,ID,same)

% MUA1: Old Mesh
% x2 , y2: nodal coordinates of new mesh
% F: scatteredInterpolant object
% Zold: Field to be mapped defined on old mesh
% Znew: Field mapped onto new mesh

Znew = zeros(length(xNew),1);

Znew(same) = F.Values(ID(same));
Znew(~same) = F(xNew(~same),yNew(~same));


end