function [UserVar,CtrlVarInRestartFile,MUA,BCs,F,l,RunInfo,Fm1]=GetInputsForForwardRestartRun(UserVar,CtrlVar,RunInfo)

narginchk(3,3) 
nargoutchk(7,8)       

fprintf('\n\n ---------  Reading restart file %s.\n',CtrlVar.NameOfRestartFiletoRead)

% For some reason the 'whos' statement sometimes fails even if the file does exist and can be loaded. So this approach, which was
% intended to make things more robust, just causes issues. 
% 
%
% try
%     Contents=whos('-file',CtrlVar.NameOfRestartFiletoRead) ;
% 
% catch exception
%     fprintf(CtrlVar.fidlog,'%s \n',exception.message);
%     error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
% end


%if any(arrayfun(@(x) isequal(x.name,'F'),Contents))
    
    try
        
        load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','RunInfo','F','l');
        
        MUAold=MUA;
        MUA=UpdateMUA(CtrlVar,MUA);
    catch exception
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
    end

    % (9 Oct 2026) Fm1, the rates of the previous time step, is saved in restart files from 9 Oct 2026 onwards. For older restart
    % files Fm1 is created with all rates set to NaN, and the explicit estimate then falls back to linear, or no, extrapolation in
    % the first time step(s) after the restart (see ExplicitEstimationUsingBackwardDifferences.m).
    WarningState=warning('off','MATLAB:load:variableNotFound');
    Stmp=load(CtrlVar.NameOfRestartFiletoRead,'Fm1');
    warning(WarningState);
    if isfield(Stmp,'Fm1') && isa(Stmp.Fm1,'UaFields')
        Fm1=Stmp.Fm1;
    else
        fprintf('\n Restart file does not contain Fm1 (an older restart file). Fm1 is created with all rates set to NaN.\n')
        Fm1=UaFields;
        Fm1.dhdt=NaN(MUAold.Nnodes,1); Fm1.dubdt=NaN(MUAold.Nnodes,1); Fm1.dvbdt=NaN(MUAold.Nnodes,1);
        Fm1.duddt=NaN(MUAold.Nnodes,1); Fm1.dvddt=NaN(MUAold.Nnodes,1);
        Fm1.dtRates=NaN;
    end
    clear Stmp
    
% else
% 
%     try
% 
%         load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','time','dt','s','b','S','B','h',...
%             'ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
%             'Itime','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
%             'GLdescriptors','l','alpha','g');
%         Co=[] ; mo=[] ; Ca=[] ; ma=[] ; dasdh=[] ; dabdh=[] ; uo=[] ; vo=[];
%         MUAold=MUA;
%         F=Vars2UaFields(ub,vb,ud,vd,uo,vo,s,b,h,S,B,AGlen,C,m,n,rho,rhow,Co,mo,Ca,ma,as,ab,dasdh,dabdh,dhdt,dsdt,dbdt,dubdt,dvbdt,duddt,dvddt,g,alpha);
% 
%         RunInfo=UaRunInfo;
% 
%     catch exception
%         fprintf(CtrlVar.fidlog,'%s \n',exception.message);
%         error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
%     end
% 
% end

F.time=CtrlVar.time ; F.dt=CtrlVar.dt ; 


% This is here to preserve past behavior from before the CtrlVar.StartTime field was introduced for the start time of the
% run. This will actually give the wrong start time in case of a restart run, but it is the best that can be done. This has
% no impact on the run results, but the waitbar might give false impression of remaining run time. 

if ~isfield(CtrlVar,"StartTime")
    CtrlVar.StartTime=CtrlVar.time;
end



% RunInfo=UaRunInfo;



if exist('BCs','var')==0
    fprintf(' The variable BCs not found in restart file. Reset. \n')
    BCs=BoundaryConditions;
end

if exist('l','var')==0
    fprintf(' The Lagrange variable l not found in restart file. Reset. \n')
    l=UaLagrangeVariables;
end

% (10 Oct 2026) The boundary conditions, including the active set of the thickness constraints (BCs.hPosNode), and the Lagrange
% multipliers as read from the restart file. Used further below to restore the active set and the multipliers.
BCsInRestartFile=BCs;
lInRestartFile=l;


if exist('RunInfo','var')==0
    fprintf(' The variable RunInfo not found in restart file. Created. \n')
    RunInfo=UaRunInfo;
end

if ~isobject(RunInfo)
   fprintf(' The variable RunInfo found in restart file is not an object. Presumably an old style restart file. \n')
   fprintf(' Recreating RunInfo as  UaRunInfo object. \n')
   RunInfo=UaRunInfo; 
end

if ~isfield(RunInfo,'Mapping') || isempty(RunInfo.Mapping)
    RunInfo.Mapping.nNewNodes=NaN;
    RunInfo.Mapping.nOldNodes=NaN;
    RunInfo.Mapping.nIdenticalNodes=NaN;
    RunInfo.Mapping.nNotIdenticalNodes=NaN;
    RunInfo.Mapping.nNotIdenticalNodesOutside=NaN;
    RunInfo.Mapping.nNotIdenticalNodesInside=NaN;
end

nRunInfo=numel(RunInfo.Forward.time) ; 
if nRunInfo < CtrlVarInRestartFile.CurrentRunStepNumber
    nRunInfo = CtrlVarInRestartFile.CurrentRunStepNumber+1000 ;
    RunInfo.Forward.time=NaN(nRunInfo,1); 
    RunInfo.Forward.dt=NaN(nRunInfo,1) ;
    RunInfo.Forward.uvhIterations=NaN(nRunInfo,1) ;
    RunInfo.Forward.uvhResidual=NaN(nRunInfo,1) ; 
    RunInfo.Forward.uvhBackTrackSteps=NaN(nRunInfo,1) ;
    RunInfo.Forward.uvhActiveSetIterations=NaN(nRunInfo,1) ;
    RunInfo.Forward.uvhActiveSetCyclical=NaN(nRunInfo,1) ;
    RunInfo.Forward.uvhActiveSetConstraints=NaN(nRunInfo,1) ;
    
end


if CtrlVar.ResetTime==1
    CtrlVarInRestartFile.time=CtrlVar.RestartTime;
    fprintf(CtrlVar.fidlog,' Time reset to CtrlVar.RestartTime=%-g \n',CtrlVarInRestartFile.time);
end


if CtrlVar.ResetTimeStep==1
    CtrlVarInRestartFile.dt=CtrlVar.dt;
    fprintf(CtrlVar.fidlog,' Time-step reset to CtrlVar.dt=%-g \n',CtrlVarInRestartFile.dt);
end


if CtrlVar.ResetRunStepNumber
    CtrlVarInRestartFile.CurrentRunStepNumber=0;
    fprintf(' RunStepNumber reset to 0 \n')
end

CtrlVar.time=CtrlVarInRestartFile.time;
CtrlVar.RestartTime=CtrlVarInRestartFile.time;
CtrlVar.dt=CtrlVarInRestartFile.dt;
CtrlVar.CurrentRunStepNumber=CtrlVarInRestartFile.CurrentRunStepNumber;

F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 

fprintf(CtrlVar.fidlog,' Starting restart run at t=%-g with dt=%-g \n',...
    CtrlVarInRestartFile.time,CtrlVarInRestartFile.dt);

if  CtrlVarInRestartFile.time> CtrlVar.EndTime
    fprintf(CtrlVar.fidlog,' Time at restart (%-g) larger than total run time (%-g) and run  is terminated. \n',CtrlVarInRestartFile.time,CtrlVar.EndTime) ;
    return
end

if CtrlVar.ReadInitialMesh==1
    fprintf(CtrlVar.fidlog,' On restart loading an initial mesh from %s \n ',CtrlVar.ReadInitialMeshFileName);
    fprintf(CtrlVar.fidlog,' This new mesh will replace the mesh in restart file. \n');
    
    
    clearvars MUA
    
    Temp=load(CtrlVar.ReadInitialMeshFileName);
    
    if isfield(Temp,'MUA')
        MUA=Temp.MUA;
        MUA=UpdateMUA(CtrlVar,MUA);
    elseif isfield(Temp,'coordinates') &&  isfield(Temp,'connectivity')
        MUA=CreateMUA(CtrlVar,Temp.connectivity,Temp.coordinates);
    else
        fprintf('Neither MUA  or connectivity and coordinates found in %s \n',CtrlVar.ReadInitialMeshFileName)
        error('Input file does not contain expected variables')
    end
    clear Temp
    
end


for I=1:CtrlVar.RefineMeshOnRestart
    fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
    
    [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    
end


isMeshChanged=HasMeshChanged(MUA,MUAold);


if isMeshChanged
    
    fprintf(CtrlVar.fidlog,' Grid changed, all variables mapped from old to new grid \n ');
    
    
    [UserVar,RunInfo,F,BCs,l]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUA,F,BCs,l);
    [RunInfo,Fm1]=MapFm1BetweenMeshes(CtrlVar,RunInfo,MUAold,MUA,Fm1);   % (9 Oct 2026)
    %[UserVar,RunInfo,F,BCs,GF]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUA,F,BCs,GF);
    
    
else
    
    if CtrlVar.TimeDependentRun
        
        % if time dependent then surface (s) and bed (b) are defined by mapping old thickness onto
        % [UserVar,~,~,F.S,F.B,F.alpha]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'SB');
        [UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,'-S-B-');
        
        l=UaLagrangeVariables;
        
    else
        
        % if a diagnostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
        fprintf('Note that as this is not a time-dependent run the ice upper and lower surfaces (s and b) are defined by the user. \n')
        fprintf('When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
        fprintf('are therefore obtained through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')
        
        
        [UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,'-s-b-S-B-rho-rhow-g');
        TestVariablesReturnedByDefineGeometryForErrors(MUA,F.s,F.b,F.S,F.B);
        %F.h=F.s-F.b;
        
    end
    
    
end

fprintf(' Note: Even though this is a restart run the following variables are defined at the beginning of the run\n')
fprintf('       through calls to corresponding user-input files: rho, rhow, g, C, m, AGlen, n, as, and ab.\n')
fprintf('       These will overwrite those in restart file.\n')


[UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);

BCs=BoundaryConditions;
[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F);

%% (10 Oct 2026) Restoring the active set and the Lagrange multipliers from the restart file
% so that a restart run continues where the previous run ended. Previously the active set was not restored (the boundary conditions are
% redefined above through DefineBoundaryConditions.m, which does not include the active set), and the multipliers were reset. The active set
% was then recovered in the first time step from h<=ThickMin (see ActiveSetInitialisation.m), but the first uvh solve started with zero
% multipliers. This is only done if the mesh has not changed (if it has, the active set is recovered by ActiveSetInitialisation.m from the
% mapped thickness) and for time-dependent runs.
%
% - The active set is restored, except for any nodes that are now fixed or tied through DefineBoundaryConditions.m.
% - The multipliers are only restored if all user-defined constraints are unchanged and the active set has been restored completely. The
%   order of the constraints, and hence of the multipliers, is then the same as in the previous run (see BCs2MLC.m).
if CtrlVar.TimeDependentRun && ~isMeshChanged

    isActiveSetRestored=true;
    nPosInFile=numel(BCsInRestartFile.hPosNode);
    if CtrlVar.ThicknessConstraints && nPosInFile>0
        [PosNode,iKeep]=setdiff(BCsInRestartFile.hPosNode(:),[BCs.hFixedNode(:);BCs.hTiedNodeA(:);BCs.hTiedNodeB(:)],'stable');
        BCs.hPosNode=PosNode;
        if numel(BCsInRestartFile.hPosValue)==nPosInFile
            BCs.hPosValue=BCsInRestartFile.hPosValue(iKeep);
            BCs.hPosValue=BCs.hPosValue(:);
        else
            BCs.hPosValue=zeros(numel(PosNode),1)+CtrlVar.ThickMin;
        end
        isActiveSetRestored=numel(PosNode)==nPosInFile;
        fprintf(' Active set restored from restart file: %i thickness constraints',numel(PosNode))
        if ~isActiveSetRestored
            fprintf(' (%i constraints in the restart file are not restored as these nodes are now fixed or tied)',nPosInFile-numel(PosNode))
        end
        fprintf('. \n')
    end

    isSame=@(a,b) isequal(a(:),b(:));
    UserConstraintsUnchanged= ...
        isSame(BCs.hFixedNode,BCsInRestartFile.hFixedNode) && isSame(BCs.hTiedNodeA,BCsInRestartFile.hTiedNodeA) && isSame(BCs.hTiedNodeB,BCsInRestartFile.hTiedNodeB) && ...
        isSame(BCs.ubFixedNode,BCsInRestartFile.ubFixedNode) && isSame(BCs.vbFixedNode,BCsInRestartFile.vbFixedNode) && ...
        isSame(BCs.ubTiedNodeA,BCsInRestartFile.ubTiedNodeA) && isSame(BCs.ubTiedNodeB,BCsInRestartFile.ubTiedNodeB) && ...
        isSame(BCs.vbTiedNodeA,BCsInRestartFile.vbTiedNodeA) && isSame(BCs.vbTiedNodeB,BCsInRestartFile.vbTiedNodeB) && ...
        isSame(BCs.ubvbFixedNormalNode,BCsInRestartFile.ubvbFixedNormalNode) && ...
        isSame(BCs.udFixedNode,BCsInRestartFile.udFixedNode) && isSame(BCs.vdFixedNode,BCsInRestartFile.vdFixedNode) && ...
        isSame(BCs.udTiedNodeA,BCsInRestartFile.udTiedNodeA) && isSame(BCs.udTiedNodeB,BCsInRestartFile.udTiedNodeB) && ...
        isSame(BCs.vdTiedNodeA,BCsInRestartFile.vdTiedNodeA) && isSame(BCs.vdTiedNodeB,BCsInRestartFile.vdTiedNodeB) && ...
        isSame(BCs.udvdFixedNormalNode,BCsInRestartFile.udvdFixedNormalNode) ;

    if UserConstraintsUnchanged && isActiveSetRestored && isa(lInRestartFile,'UaLagrangeVariables')
        l=lInRestartFile;
        fprintf(' Lagrange multipliers restored from restart file. \n')
    elseif ~UserConstraintsUnchanged
        fprintf(' Lagrange multipliers not restored from restart file, because the boundary conditions defined in DefineBoundaryConditions.m have changed. \n')
    elseif ~isActiveSetRestored
        fprintf(' Lagrange multipliers not restored from restart file, because the active set could not be fully restored. \n')
    end

end
% (10 Oct 2026) Optional reset of the accumulated time-discretisation error estimate at a restart (by default it is continued).
if isfield(CtrlVar,"TimeDiscretisationErrorEstimate") && isfield(CtrlVar.TimeDiscretisationErrorEstimate,"ResetAccumulatedAtRestart") ...
        && CtrlVar.TimeDiscretisationErrorEstimate.ResetAccumulatedAtRestart && ~isempty(F.hTimeDiscretisationErrorAccumulated)
    F.hTimeDiscretisationErrorAccumulated=[];
    fprintf(' Accumulated time-discretisation error estimate reset at restart. \n')
end

% This is now a part of DefineSlipperiness
% if CtrlVar.IncludeMelangeModelPhysics
%     fprintf(' Also here defining Melange/Sea-ice model parameters through a call to a user-input file. \n')
%     [UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,F);
% end
% 

if CtrlVar.doplots==1 && CtrlVar.PlotBCs==1
    
    fig=FindOrCreateFigure("Boundary Conditions");
    clf(fig) 
    hold off
    PlotBoundaryConditions(CtrlVar,MUA,BCs);
    
end




fprintf(' ---------   Reading restart file and defining start values for restart run is now done.\n\n')




end