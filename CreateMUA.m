function MUA=CreateMUA(CtrlVar,connectivity,coordinates,CalcMUA_Derivatives,FindMUA_Boundary)

% MUA=CreateMUA(CtrlVar,connectivity,coordinates,CalcMUA_Derivatives,FindMUA_Boundary)
%
% Creates the Úa mesh structure containing all information about the FE mesh
% such as coordinates, connectivity, boundary nodes, etc Also (optionally)
% calculates element derivatives used in the matrix assembly.
%
% Example: MUA=CreateMUA(CtrlVar,connectivity,coordinates);
%
%

if nargin<4
    CalcMUA_Derivatives=1;
end

if nargin<5;
    FindMUA_Boundary=1;
end

% eliminate coordinates that are not part of mesh, and update connectivity accordingly
[K,~,J]=unique(connectivity(:));
connectivity=reshape(J,size(connectivity));
coordinates=coordinates(K,:);

% First check if element type on input is as reqested by user, and if not change
[MUA.coordinates,MUA.connectivity]=ChangeElementType(coordinates,connectivity,CtrlVar.TriNodes);
MUA.Nnodes=size(MUA.coordinates,1);
MUA.Nele=size(MUA.connectivity,1);
MUA.nod=size(MUA.connectivity,2);

CtrlVar=NrOfIntegrationPoints(CtrlVar);
MUA.nip=CtrlVar.nip ;
MUA.niph=CtrlVar.niph;


ndim=2;
[MUA.points,MUA.weights]=sample('triangle',MUA.nip,ndim);


if FindMUA_Boundary
    [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
end

if CalcMUA_Derivatives
    [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates);
end

end