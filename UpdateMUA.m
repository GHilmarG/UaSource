function MUA=UpdateMUA(CtrlVar,MUA)

%%
% MUA=UpdateMUA(CtrlVar,MUA)
% Updates MUA and calculates any missing fields.
% On input MUA must have the fields coordinates and connectivity.
%

if ~isfield(MUA,'coordinates')
    error('MUA must have a coordinates field')
end

if ~isfield(MUA,'connectivity')
    error('MUA must have a connectivity field')
end


% first make sure that the element is of the right type
[MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,CtrlVar.TriNodes);

% checking if MUA has the QuadratureRuleDegree field
if CtrlVar.QuadRules2021  && ~isfield(MUA,"QuadratureRuleDegree")
    QuadratureFieldMissing=true;
else
    QuadratureFieldMissing=false ;
end



if ~isfield(MUA,'niph')  || ~isfield(MUA,'nip')  ||  ~isfield(MUA,'points')  || ~isfield(MUA,'weights') ||  QuadratureFieldMissing 

    if CtrlVar.QuadRules2021
        Degree=QuadratureRuleDegree(CtrlVar);
        Q=quadtriangle(Degree,'Type','nonproduct','Points','inside','Domain',[0 0 ; 1 0 ; 0 1]) ;
        
        MUA.QuadratureRuleDegree=Degree;
        MUA.nip=size(Q.Points,1);
        MUA.niph=size(Q.Points,1);
        MUA.points=Q.Points;
        MUA.weights=Q.Weights;
        
    else
        CtrlVar=NrOfIntegrationPoints(CtrlVar);
        
        MUA.QuadratureRuleDegree=nan;
        MUA.niph=CtrlVar.niph;
        MUA.nip=CtrlVar.nip;
        [MUA.points,MUA.weights]=sample('triangle',MUA.nip,MUA.ndim);
    end

end

if ~isfield(CtrlVar,'MUA')
    CtrlVar.MUA.MassMatrix=0;
    CtrlVar.MUA.StiffnessMatrix=0;
end


if ~isfield(CtrlVar.MUA,'MassMatrix')
    CtrlVar.MUA.MassMatrix=0;
end


if ~isfield(CtrlVar.MUA,'StiffnessMatrix')
    CtrlVar.MUA.StiffnessMatrix=0;
end

if  ~isfield(CtrlVar,'InfoLevel')
    CtrlVar.InfoLevel=1;
end

if ~isfield(CtrlVar,'FindMUA_Boundary')
    CtrlVar.FindMUA_Boundary=0;
end

if ~isfield(CtrlVar,'CalcMUA_Derivatives')
    CtrlVar.CalcMUA_Derivatives=0;
end


if CtrlVar.FindMUA_Boundary && isempty(MUA.TR)
    [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
end


%% Now consider the possibility that we are using the post 2021 quad rules and that the quadrature degree has changed

 QuadratureRuleHasChanged=false;

if CtrlVar.QuadRules2021
    Degree=QuadratureRuleDegree(CtrlVar);
    if MUA.QuadratureRuleDegree~=Degree
        QuadratureRuleHasChanged=true;
    end
end

%% Now consider the possibility the FE coordinates and connectivity has changed
% and that the other fields are not up to date

MeshHasChanged = ...
    MUA.nod~=size(MUA.connectivity,2) || ...
    MUA.Nele~=size(MUA.connectivity,1) || ...
    MUA.Nnodes~=size(MUA.coordinates,1) || ...
    QuadratureRuleHasChanged ; 

if MeshHasChanged
    if CtrlVar.InfoLevel>=10
        fprintf('UpdateMUA: Mesh has changed \n ')
        fprintf('UpdateMUA: finding mesh boundary \n ')
        fprintf('UpdateMUA: Calculating mesh derivatives \n ')
    end
    MUA.ndim=2;
    MUA.nod=size(MUA.connectivity,2);
    MUA.Nele=size(MUA.connectivity,1);
    MUA.Nnodes=size(MUA.coordinates,1);
    
    if CtrlVar.QuadRules2021
        
        Degree=QuadratureRuleDegree(CtrlVar);
        MUA.QuadratureRuleDegree=Degree;
        Q=quadtriangle(Degree,'Type','nonproduct','Points','inside','Domain',[0 0 ; 1 0 ; 0 1]) ;
        MUA.nip=size(Q.Points,1);
        MUA.niph=size(Q.Points,1);
        MUA.points=Q.Points;
        MUA.weights=Q.Weights;

    else
        
        CtrlVar=NrOfIntegrationPoints(CtrlVar);
        MUA.QuadratureRuleDegree=nan;
        MUA.nip=CtrlVar.nip ; MUA.niph=CtrlVar.niph;
        [MUA.points,MUA.weights]=sample('triangle',MUA.nip,MUA.ndim);
    end


    if CtrlVar.FindMUA_Boundary
        [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
    else
        MUA.Boundary=[];
        MUA.TR=[];
    end
    
    if CtrlVar.CalcMUA_Derivatives
        [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
    else
        MUA.Deriv=[];
        MUA.DetJ=[];
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
    %     MUA.L=chol(MUA.M,'upper');
    % end


    [MUA.xEle,MUA.yEle]=ElementCoordinates(MUA.connectivity,MUA.coordinates);

    MUA.workers=[];

    if CtrlVar.Parallel.uvAssembly.spmd.isOn || CtrlVar.Parallel.uvhAssembly.spmd.isOn

        MUA.workers=BuildMuaWorkers(CtrlVar,MUA,MUA.workers) ;

    end


end



%% and now the possibility that the mesh has not changed but some of the fields were not defined previously
if  CtrlVar.FindMUA_Boundary
    if ~isfield(MUA,'Boundary') || ~isfield(MUA.Boundary,'x') || ~isfield(MUA.Boundary,'y')
        fprintf('UpdateMUA: finding mesh bounday \n ')
        MUA.Boundary=FindBoundary(MUA.connectivity,MUA.coordinates);
    end
end



if CtrlVar.CalcMUA_Derivatives
    if ~isfield(MUA,'DetJ') || ~isfield(MUA,'Deriv')
        
        [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
    end
end

[NeleTest,ndimTest,nodTest,nipTest]=size(MUA.Deriv);
MUADerivHasChanged=isempty(MUA.Deriv)  ||  NeleTest~=MUA.Nele || nodTest~=MUA.nod || nipTest~=MUA.nip ; 


if CtrlVar.CalcMUA_Derivatives && MUADerivHasChanged
        [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end


if  (CtrlVar.MUA.MassMatrix || CtrlVar.MUA.DecomposeMassMatrix ) &&  ( ~isfield(MUA,'M') || isempty(MUA.M) || MUADerivHasChanged  )  
    MUA.M=MassMatrix2D1dof(MUA);
end


if CtrlVar.MUA.DecomposeMassMatrix  && ( ~isfield(MUA,'dM')  ||isempty(MUA.dM)  || ~all(MUA.dM.MatrixSize==size(MUA.M)) || MUADerivHasChanged)
    MUA.dM=decomposition(MUA.M,'chol','upper') ;
end


if CtrlVar.MUA.StiffnessMatrix &&  (~isfield(MUA,'Dxx')  || MUADerivHasChanged)
    [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
end

% if CtrlVar.Inverse.AdjointGradientPreMultiplier=="M"
%    MUA.L=chol(MUA.M,'upper');
% end

if ~isfield(MUA,'xEle')
    [MUA.xEle,MUA.yEle]=ElementCoordinates(MUA.connectivity,MUA.coordinates);
end


if ~isfield(MUA,'Boundary') ||  ~isfield(MUA,'TR')
    if CtrlVar.FindMUA_Boundary
        [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
    else
        MUA.Boundary=[];
        MUA.TR=[];
    end
end

MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity); % areas if each element
MUA.Area=sum(MUA.EleAreas);                               % total FE mesh area


if ~isfield(MUA,"workers")
    MUA.workers=[];
end



if ( CtrlVar.Parallel.uvAssembly.spmd.isOn || CtrlVar.Parallel.uvhAssembly.spmd.isOn  )

    % poolobj = gcp;
    poolobj = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already

    if isempty(poolobj)

        fprintf("SPMD assembly is set to true, but parallel pool is empty. \n")
        fprintf(" Create a parallel pool ahead of the call to %ca.\n",218)

    else
        CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;

        % not sure how to best to check if the composite is in good state
        % if length(MUA.workers)  ~= CtrlVar.Parallel.uvhAssembly.spmd.nWorkers
        %     MUA.workers=[] ;
        % end

        if ~all(exist(MUA.workers))
            MUA.workers=[];
        end

        MUA.workers=BuildMuaWorkers(CtrlVar,MUA,MUA.workers) ;
    end

end



end