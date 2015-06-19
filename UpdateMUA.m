function MUA=UpdateMUA(CtrlVar,MUA)

%%
% MUA=UpdateMUA(CtrlVar,MUA)
% This updates MUA and calculates any missing fields
% On input MUA must have the fields coordinates and connectivity

if ~isfield(MUA,'coordinates')
    error('MUA must have a coordinates field')
end

if ~isfield(MUA,'connectivity')
    error('MUA must have a connectivity field')
end

% first make sure that the element is of the right type
[MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,CtrlVar.TriNodes);


CtrlVar=NrOfIntegrationPoints(CtrlVar);


%% Now consider the possibility the FE coordinates and connectivity has changed
% and that the other fields are not up to date

MeshHasChanged = ...
    MUA.nod~=size(MUA.connectivity,2) || ...
    MUA.Nele~=size(MUA.connectivity,1) || ...
    MUA.Nnodes~=size(MUA.coordinates,1) || ...
    MUA.nip~=CtrlVar.nip && MUA.niph~=CtrlVar.niph ;

%    all(MUA.points==points) && ...
%    all(MUA.weights==weights) ;Ct

if MeshHasChanged
    if CtrlVar.InfoLevel>=10;
        fprintf('UpdateMUA: Mesh has changed \n ')
        fprintf('UpdateMUA: finding mesh boundary \n ')
        fprintf('UpdateMUA: Calculating mesh derivatives \n ')
    end
    MUA.ndim=2;
    
    CtrlVar=NrOfIntegrationPoints(CtrlVar); 
    MUA.nip=CtrlVar.nip ; MUA.niph=CtrlVar.niph;
    
    MUA.Nele=size(MUA.connectivity,1);
    MUA.Nnodes=size(MUA.coordinates,1);
    MUA.nod=size(MUA.connectivity,2);
    [MUA.points,MUA.weights]=sample('triangle',MUA.nip,MUA.ndim);
    
    [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
    [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates);
    
end



%% and now the possibility that the mesh has not changed but some of the fields were not defined previously
if ~isfield(MUA,'Boundary')
    fprintf('UpdateMUA: finding mesh bounday \n ')
    
    MUA.Boundary=FindBoundary(MUA.connectivity,MUA.coordinates);
end

if ~isfield(MUA,'points')  || ~isfield(MUA,'weights')
    MUA.ndim=2;
    CtrlVar=NrOfIntegrationPoints(CtrlVar); MUA.nip=CtrlVar.nip ; MUA.niph=CtrlVar.niph;
    [MUA.points,MUA.weights]=sample('triangle',MUA.nip,MUA.ndim);
end

if ~isfield(MUA,'DetJ') || ~isfield(MUA,'Deriv')
    [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates);
end


[NeleTest,ndimTest,nodTest,nipTest]=size(MUA.Deriv);

if NeleTest~=MUA.Nele || nodTest~=MUA.nod || nipTest~=MUA.nip
    [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates);
end


end