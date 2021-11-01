function [UserVar,CtrlVarInRestartFile,MUA,BCs,F,l,RunInfo]=GetInputsForForwardRestartRun(UserVar,CtrlVar,RunInfo)

narginchk(3,3) 
nargoutchk(7,7)       

fprintf('\n\n ---------  Reading restart file %s.\n',CtrlVar.NameOfRestartFiletoRead)

try
    Contents=whos('-file',CtrlVar.NameOfRestartFiletoRead) ;
    
catch exception
    fprintf(CtrlVar.fidlog,'%s \n',exception.message);
    error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
end


if any(arrayfun(@(x) isequal(x.name,'F'),Contents))
    
    try
        
        load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','RunInfo','F','l');
        
        MUAold=MUA;
        MUA=UpdateMUA(CtrlVar,MUA);
    catch exception
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
    end
    
else
    
    try
        
        load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','time','dt','s','b','S','B','h',...
            'ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
            'Itime','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
            'GLdescriptors','l','alpha','g');
        Co=[] ; mo=[] ; Ca=[] ; ma=[] ; dasdh=[] ; dabdh=[] ; uo=[] ; vo=[];
        MUAold=MUA;
        F=Vars2UaFields(ub,vb,ud,vd,uo,vo,s,b,h,S,B,AGlen,C,m,n,rho,rhow,Co,mo,Ca,ma,as,ab,dasdh,dabdh,dhdt,dsdt,dbdt,dubdt,dvbdt,duddt,dvddt,g,alpha);
        
        RunInfo=UaRunInfo;
        
    catch exception
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
    end
    
end

F.time=CtrlVar.time ; F.dt=CtrlVar.dt ; 


% RunInfo=UaRunInfo;
RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";
if CtrlVar.WriteRunInfoFile
    RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
end

if exist('BCs','var')==0
    fprintf(' The variable BCs not found in restart file. Reset. \n')
    BCs=BoundaryConditions;
end

if exist('l','var')==0
    fprintf(' The Lagrange variable l not found in restart file. Reset. \n')
    l=UaLagrangeVariables;
end


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

if  CtrlVarInRestartFile.time> CtrlVar.TotalTime
    fprintf(CtrlVar.fidlog,' Time at restart (%-g) larger than total run time (%-g) and run  is terminated. \n',CtrlVarInRestartFile.time,CtrlVar.TotalTime) ;
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


%[UserVar,F]=GetDensities(UserVar,CtrlVar,MUA,F);
%[F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
[UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);

BCs=BoundaryConditions;
[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F);

if CtrlVar.IncludeMelangeModelPhysics
    fprintf(' Also here defining Melange/Sea-ice model parameters through a call to a user-input file. \n')
    [UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,BCs,F);
end


if CtrlVar.doplots==1 && CtrlVar.PlotBCs==1
    
    fig=FindOrCreateFigure("Boundary Conditions");
    clf(fig) 
    hold off
    PlotBoundaryConditions(CtrlVar,MUA,BCs);
    
end




fprintf(' ---------   Reading restart file and defining start values for restart run is now done.\n\n')




end