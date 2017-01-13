function [UserVar,CtrlVarInRestartFile,MUA,BCs,F,l,RunInfo]=GetInputsForForwardRestartRun(UserVar,CtrlVar)


fprintf('\n\n ---------  Reading restart file and defining start values for restart run.\n')

RunInfo=[];

Contents=whos('-file',CtrlVar.NameOfRestartFiletoRead) ;

if any(arrayfun(@(x) isequal(x.name,'F'),Contents))
    
    try
        
        load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','RunInfo','time','dt','F','GF','l');
        
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
        
        
        F=Vars2UaFields(ub,vb,ud,vd,uo,vo,s,b,h,S,B,AGlen,C,m,n,rho,rhow,Co,mo,Ca,ma,as,ab,dasdh,dabdh,dhdt,dsdt,dbdt,dubdt,dvbdt,duddt,dvddt,g,alpha);
        
    catch exception
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
    end
    
end


% Thickness should only depend on s and b in restart file
% (The only exeption being that if h is less than CtrlVar.ThickMin,
% and CtrlVar.ResetThicknessToMinThickness true, then h is first modified accordingly.)
F.h=F.s-F.b;

[F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
%[F.b,F.s,F.h]=Calc_bs_From_hBS(F.h,F.S,F.B,F.rho,F.rhow,CtrlVar,MUA.coordinates);

if exist('MUA','var')==0
    fprintf(' The variable MUA not found in restart file. Try to read connectivity and coordinates from restart file and then to create MUA \n')
    load(CtrlVar.NameOfRestartFiletoRead,'connectivity','coordinates')
    MUA=CreateMUA(CtrlVar,connectivity,coordinates,1,1);
end

if exist('BCs','var')==0
    fprintf(' The variable BCs not found in restart file. Reset. \n')
    BCs=BoundaryConditions;
end

if exist('l','var')==0
    fprintf(' The Lagrange variable l not found in restart file. Reset. \n')
    l=UaLagrangeVariables;
end


if CtrlVar.ResetTime==1
    CtrlVarInRestartFile.time=CtrlVar.time;
    CtrlVarInRestartFile.CurrentRunStepNumber=0;
    fprintf(CtrlVar.fidlog,' Time reset to %-g \n',time);
end

if CtrlVar.ResetTimeStep==1
    CtrlVarInRestartFile.dt=CtrlVar.dt;
    fprintf(CtrlVar.fidlog,' Time-step reset to %-g \n',dt);
end

fprintf(CtrlVar.fidlog,' Read restart file %s.  Starting restart run at t=%-g with dt=%-g \n',CtrlVar.NameOfRestartFiletoRead,time,dt);

if  CtrlVarInRestartFile.time> CtrlVar.TotalTime
    fprintf(CtrlVar.fidlog,' Time at restart (%-g) larger than total run time (%-g) and run  is terminated. \n',time,CtrlVar.TotalTime) ;
    return
end

MUAold=MUA;


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

MeshChanged=HasMeshChanged(MUAold,MUA);

if MeshChanged
    
    fprintf(CtrlVar.fidlog,' Grid changed, all variables mapped from old to new grid \n ');
    
    [UserVar,F,BCs,GF]=MapFbetweenMeshes(UserVar,CtrlVar,MUAold,MUA,F,BCs,GF);
    
    
else
    
    if CtrlVar.TimeDependentRun
        
        % if time dependent then surface (s) and bed (b) are defined by mapping old thickness onto
        [UserVar,~,~,F.S,F.B,F.alpha]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'SB');
        
        l=UaLagrangeVariables;
        
    else
        
        % if a diagnostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
        fprintf('Note that as this is not a time-dependent run the ice upper and lower surfaces (s and b) are defined by the user. \n')
        fprintf('When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
        fprintf('are therefore obtained through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')
        
        [UserVar,F.s,F.b,F.S,F.B,F.alpha]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'sbSB');
        TestVariablesReturnedByDefineGeometryForErrors(MUA,F.s,F.b,F.S,F.B);
        F.h=F.s-F.b;
        
    end
    
    
end

fprintf(' Note: Even though this is a restart run the following variables are defined at the beginning of the run\n')
fprintf('       through calls to corresponding user-input files: rho, rhow, g, C, m, AGlen, n, as, and ab.\n')
fprintf('       These will owerwrite those in restart file.\n')


[UserVar,F]=GetDensities(UserVar,CtrlVar,MUA,F);
[F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);

%[F.b,F.s,F.h]=Calc_bs_From_hBS(F.h,F.S,F.B,F.rho,F.rhow,CtrlVar,MUA.coordinates);
%GF = GL2d(F.B,F.S,F.h,F.rhow,F.rho,MUA.connectivity,CtrlVar);

[UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F,GF);
[UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F,GF);
[UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);

BCs=BoundaryConditions;
[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F,GF);

if CtrlVar.IncludeMelangeModelPhysics
    fprintf(' Also here defining Melange/Sea-ice model parameters through a call to a user-input file. \n')
    [UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,BCs,F,GF);
end


if CtrlVar.doplots==1 && CtrlVar.PlotBCs==1
    
    figure
    PlotBoundaryConditions(CtrlVar,MUA,BCs);
    
end




fprintf(' ---------   Reading restart file and defining start values for restart run is now done.\n\n')




end