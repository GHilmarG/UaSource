

function MUA=UpdateMUA(CtrlVar,MUA)

%%
% MUA=UpdateMUA(CtrlVar,MUA)
%
% Updates MUA and calculates any missing fields.
%
% On input MUA must have the fields coordinates and connectivity.
%
% This is typically called whenever the FE-mesh changes. It
%

if ~isfield(MUA,'coordinates')
    error('MUA must have a coordinates field')
end

if ~isfield(MUA,'connectivity')
    error('MUA must have a connectivity field')
end

%% Start doing some simple overall checks related to quadrature rules which were changed in 2021.
% This is here for mainly completeness and so that users can re-run old examples and use old re-start files.

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

QuadratureRuleHasChanged=false;
% has the quadrature degree changed?
if CtrlVar.QuadRules2021 &&  ~isempty(CtrlVar.QuadratureRuleDegree)  &&  MUA.QuadratureRuleDegree ~= CtrlVar.QuadratureRuleDegree

    Degree=QuadratureRuleDegree(CtrlVar);
    Q=quadtriangle(Degree,'Type','nonproduct','Points','inside','Domain',[0 0 ; 1 0 ; 0 1]) ;

    MUA.QuadratureRuleDegree=Degree;
    MUA.nip=size(Q.Points,1);
    MUA.niph=size(Q.Points,1);
    MUA.points=Q.Points;
    MUA.weights=Q.Weights;
    QuadratureRuleHasChanged=true;

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



if CtrlVar.QuadRules2021
    Degree=QuadratureRuleDegree(CtrlVar);
    if MUA.QuadratureRuleDegree~=Degree
        QuadratureRuleHasChanged=true;
    end
end

%% Now check if the mesh has changed.
% (17 Oct 2026) In addition to the numbers of nodes and elements, the key of the mesh (see MeshKey.m) is compared with the one stored
% in MUA.MeshKey when the quantities derived from the mesh were last calculated. TopologyHasChanged: the connectivity has changed
% (including a renumbering of nodes or elements). GeometryHasChanged: the coordinates have changed. If MUA has no key (eg MUA from an
% older file), everything is recalculated once.

NewMeshKey=MeshKey(MUA) ;
if isfield(MUA,"MeshKey") && isstruct(MUA.MeshKey) && isfield(MUA.MeshKey,"Topology") && isfield(MUA.MeshKey,"Geometry")
    TopologyHasChanged=~isequal(MUA.MeshKey.Topology,NewMeshKey.Topology) ;
    GeometryHasChanged=~isequal(MUA.MeshKey.Geometry,NewMeshKey.Geometry) ;
else
    TopologyHasChanged=true ;
    GeometryHasChanged=true ;
end

SizesHaveChanged = ...
    MUA.nod~=size(MUA.connectivity,2) || ...
    MUA.Nele~=size(MUA.connectivity,1) || ...
    MUA.Nnodes~=size(MUA.coordinates,1) ;

TopologyHasChanged = TopologyHasChanged || SizesHaveChanged ;
MeshHasChanged = TopologyHasChanged || GeometryHasChanged || QuadratureRuleHasChanged ;


% If the mesh has changed, quite a few things may have to be re-calculated. This is done here below, and this is the main
% objective of this function.
if MeshHasChanged
    if CtrlVar.InfoLevel>=10
        fprintf('UpdateMUA: Mesh has changed \n ')
        fprintf('UpdateMUA: topology changed %d, geometry changed %d, quadrature rule changed %d \n ',TopologyHasChanged,GeometryHasChanged,QuadratureRuleHasChanged)
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


    if ~isfield(MUA,"uvAssemblyPattern")
        MUA.uvAssemblyPattern=[];
    end
    if ~isfield(MUA,"uvhAssemblyPattern")
        MUA.uvhAssemblyPattern=[];
    end
    if ~isfield(MUA,"hAssemblyPattern")
        MUA.hAssemblyPattern=[];
    end

    if (CtrlVar.MUA.AssemblyPattern.uv || CtrlVar.MUA.AssemblyPattern.uvh || (isfield(CtrlVar.MUA.AssemblyPattern,"h") && CtrlVar.MUA.AssemblyPattern.h)) && TopologyHasChanged   % (17 Oct 2026) the patterns only depend on the connectivity
        [MUA.uvAssemblyPattern,MUA.uvhAssemblyPattern,MUA.hAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);   % (31 Oct 2026) and the h pattern
    end



    CholeskyMassMatrix=isfield(CtrlVar.MUA,"CholeskyMassMatrix") && CtrlVar.MUA.CholeskyMassMatrix ;
    if CtrlVar.MUA.MassMatrix || CtrlVar.MUA.DecomposeMassMatrix || CholeskyMassMatrix

        MUA.M=MassMatrix2D1dof(MUA);

        if CtrlVar.MUA.DecomposeMassMatrix
            MUA.dM=decomposition(MUA.M,'chol','upper') ;
        end

        if CholeskyMassMatrix   % (17 Oct 2026) previously not recalculated here, so that MUA.MC and MUA.Mp were left unchanged
            [MUA.MC,~,MUA.Mp]=chol(MUA.M,"vector");
        end



    else

        MUA.M=[] ; MUA.dM=[] ;

    end


    if CtrlVar.MUA.StiffnessMatrix
        [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
    end



    [MUA.xEle,MUA.yEle]=ElementCoordinates(MUA.connectivity,MUA.coordinates);

    % (17 Oct 2026) the element areas are needed by BuildMuaWorkers (and the SUPG terms of the assembly use them), so they have to be
    % recalculated before the workers are built; previously they were only recalculated further below, so the workers got the old areas
    MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity);
    MUA.Area=sum(MUA.EleAreas);

    MUA.workers=[];      % (26 Oct 2026) the worker data are built anew at the end of this function (see there)
    MUA.workersKey=[];


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

    if CtrlVar.MUA.DecomposeMassMatrix
        MUA.dM=decomposition(MUA.M,'chol','upper') ;
    end

    if isfield(CtrlVar.MUA,"CholeskyMassMatrix") && CtrlVar.MUA.CholeskyMassMatrix

        [MUA.MC,~,MUA.Mp]=chol(MUA.M,"vector");

    end


end

if ~isfield(MUA,"uvAssemblyPattern")
    MUA.uvAssemblyPattern=[];
end
if ~isfield(MUA,"uvhAssemblyPattern")
    MUA.uvhAssemblyPattern=[];
end
if ~isfield(MUA,"hAssemblyPattern")
    MUA.hAssemblyPattern=[];
end

if ~isfield(CtrlVar,"MUA")  || ~isfield(CtrlVar.MUA,"AssemblyPattern")
    CtrlVar.MUA.AssemblyPattern.uv=true;
    CtrlVar.MUA.AssemblyPattern.uvh=true;
end

if ~isfield(CtrlVar.MUA.AssemblyPattern,"uv")
    CtrlVar.MUA.AssemblyPattern.uv=true;
end

if ~isfield(CtrlVar.MUA.AssemblyPattern,"uvh")
    CtrlVar.MUA.AssemblyPattern.uvh=true;
end

if ~isfield(CtrlVar.MUA.AssemblyPattern,"h")
    CtrlVar.MUA.AssemblyPattern.h=true;
end



if CtrlVar.MUA.AssemblyPattern.uv && CtrlVar.MUA.AssemblyPattern.uvh
    if isempty(MUA.uvAssemblyPattern)  &&  isempty(MUA.uvhAssemblyPattern)
        [MUA.uvAssemblyPattern,MUA.uvhAssemblyPattern,MUA.hAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);
    end
end


if CtrlVar.MUA.AssemblyPattern.uv && isempty(MUA.uvAssemblyPattern)
    MUA.uvAssemblyPattern=AssemblyPatternCache(CtrlVar,MUA);
end


if CtrlVar.MUA.AssemblyPattern.uvh && isempty(MUA.uvhAssemblyPattern)
    [~,MUA.uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);
end

if CtrlVar.MUA.AssemblyPattern.h && isempty(MUA.hAssemblyPattern)   % (31 Oct 2026), eg a mesh from a file written before the h pattern existed
    [~,~,MUA.hAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);
end


%% it is possible that the mesh has not changed, but the AssemblyPattern is outdated.
if CtrlVar.MUA.AssemblyPattern.uv  && MUA.uvAssemblyPattern.neq ~= 2* MUA.Nnodes

 
    MUA.uvAssemblyPattern=AssemblyPatternCache(CtrlVar,MUA);

end

if CtrlVar.MUA.AssemblyPattern.uvh  && MUA.uvhAssemblyPattern.neq ~= 3* MUA.Nnodes

  
    [~,MUA.uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);

end

if CtrlVar.MUA.AssemblyPattern.h  && MUA.hAssemblyPattern.neq ~= MUA.Nnodes
    [~,~,MUA.hAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA);
end

%% It is possible that the decomposition object has somehow become invalid. Not sure how, but if, for example a mesh is re-read then possibly the decomposition object is still there but invalid
%
if CtrlVar.MUA.DecomposeMassMatrix  &&   ( ~isfield(MUA,'dM') || isempty(MUA.dM)  || any(MUA.dM.MatrixSize==[0 0]))

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

% (26 Oct 2026) the element areas are recalculated above if the mesh has changed; here only if missing
if ~isfield(MUA,"EleAreas") || numel(MUA.EleAreas)~=MUA.Nele || ~isfield(MUA,"Area")
    MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity); % areas if each element
    MUA.Area=sum(MUA.EleAreas);                               % total FE mesh area
end


if ~isfield(MUA,"workers")
    MUA.workers=[];
end
if ~isfield(MUA,"workersKey")
    MUA.workersKey=[];
end

if ( CtrlVar.Parallel.uvAssembly.spmd.isOn || CtrlVar.Parallel.uvhAssembly.spmd.isOn  )

    poolobj = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already

    if isempty(poolobj)

        fprintf("SPMD assembly is set to true, but parallel pool is empty. Create a parallel pool ahead of the call to %ca.\n",218)

    else
        CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;

        % (26 Oct 2026) Are the worker data (the Composite MUA.workers) valid for this mesh and this pool?
        %  - all(exist(MUA.workers)) is false for all workers if the pool has been deleted, or replaced by a new one (even of the same
        %    size): it checks that the Composite belongs to the current pool
        %  - numel(MUA.workers) is the number of workers of the pool for which the worker data were built
        %  - the contents of a Composite can only be checked by communicating with the workers, so the key of the mesh, the element
        %    type and quadrature rule, and the number of workers for which they were built are kept in MUA.workersKey (see
        %    MuaWorkersKey.m) and compared with the current ones. This also detects a change of the mesh that leaves the numbers of
        %    nodes and elements unchanged.
        % This takes well under a millisecond. Only if one of these does not hold are the worker data built (anew). Previously,
        % BuildMuaWorkers.m checked this itself in each call, with an spmd round trip and a transfer of the worker data to the client
        % (about 45 ms per call for a mesh of 110,000 elements).
        Key=MuaWorkersKey(NewMeshKey,MUA,poolobj.NumWorkers) ;
        WorkersValid = ~isempty(MUA.workers) && isequal(MUA.workersKey,Key) && numel(MUA.workers)==poolobj.NumWorkers && all(exist(MUA.workers)) ;  %#ok<EXIST>

        if ~WorkersValid
            MUA.workers=BuildMuaWorkers(CtrlVar,MUA,[]) ;
            if isempty(MUA.workers)        % not built (CtrlVar.Parallel.BuildWorkers=false)
                MUA.workersKey=[] ;
            else
                MUA.workersKey=Key ;
            end
        end
    end

end




% (17 Oct 2026) the key of the mesh for which the quantities derived from the mesh were calculated (see MeshKey.m)
MUA.MeshKey=NewMeshKey;

end