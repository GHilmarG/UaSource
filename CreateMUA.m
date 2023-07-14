function MUA=CreateMUA(CtrlVar,connectivity,coordinates,RefineMesh)

% MUA=CreateMUA(CtrlVar,connectivity,coordinates,CalcMUA_Derivatives,FindMUA_Boundary)
%
% Creates the Úa mesh structure containing all information about the FE mesh
% such as coordinates, connectivity, boundary nodes, etc Also (optionally)
% calculates element derivatives used in the matrix assembly.
%
% Example: MUA=CreateMUA(CtrlVar,connectivity,coordinates);
%
%

narginchk(3,4)

if nargin<4
    RefineMesh=[];
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

if CtrlVar.QuadRules2021
    % This uses the new quad rules implemented in 2021
    % The old ones are there for comparision and compatability
    Degree=QuadratureRuleDegree(CtrlVar);
    MUA.QuadratureRuleDegree=Degree;
    Q=quadtriangle(Degree,'Type','nonproduct','Points','inside','Domain',[0 0 ; 1 0 ; 0 1]) ;
    MUA.nip=size(Q.Points,1);
    MUA.niph=size(Q.Points,1);
    CtrlVar.nip=MUA.nip;
    CtrlVar.niph=MUA.niph;
    MUA.points=Q.Points;
    MUA.weights=Q.Weights;
else
    CtrlVar=NrOfIntegrationPoints(CtrlVar);
    MUA.QuadratureRuleDegree=nan;
    MUA.nip=CtrlVar.nip ;
    MUA.niph=CtrlVar.niph;
    ndim=2;
    [MUA.points,MUA.weights]=sample('triangle',MUA.nip,ndim);
end


[MUA.connectivity]=TestAndCorrectForInsideOutElements(CtrlVar,MUA.coordinates,MUA.connectivity);

if CtrlVar.CalcMUA_Derivatives
    [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates,MUA.nip,MUA.points);
else
    MUA.Deriv=[];
    MUA.DetJ=[];
end

if CtrlVar.FindMUA_Boundary
    [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
else
    MUA.Boundary=[];
    MUA.TR=[];
end

if CtrlVar.MUA.MassMatrix || CtrlVar.MUA.DecomposeMassMatrix
    MUA.M=MassMatrix2D1dof(MUA);
end

if CtrlVar.MUA.DecomposeMassMatrix
    MUA.dM=decomposition(MUA.M,'chol','upper') ;  
end


if CtrlVar.MUA.StiffnessMatrix
    [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
end

% if CtrlVar.Inverse.AdjointGradientPreMultiplier=="M"
%    MUA.L=chol(MUA.M,'upper');
% end

[MUA.xEle,MUA.yEle]=ElementCoordinates(MUA.connectivity,MUA.coordinates);

if ~isempty(RefineMesh)
    MUA.RefineMesh=RefineMesh;
end

MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity); % areas for each element
MUA.Area=sum(MUA.EleAreas);                               % total FE mesh area





end