function    [MUA,BCs,time,dt,CurrentRunStepNumber,s,b,S,B,ub,vb,ud,vd,l,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
    dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
    GF,GLdescriptors]=GetInputsForForwardRestartRun(CtrlVar)


try
    
    load(CtrlVar.NameOfRestartFiletoRead,'CtrlVarInRestartFile','MUA','BCs','time','dt','s','b','S','B','h','ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
        'Itime','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
        'GLdescriptors','l');
    
catch exception
    fprintf(CtrlVar.fidlog,'%s \n',exception.message);
    error('could not load restart file %s',CtrlVar.NameOfRestartFiletoRead)
end

CurrentRunStepNumber=Itime;  % I used to refer to this as Itime, change later

% Thickness should only depend on s and b in restart file
% (The only exeption being that if h is less than CtrlVar.ThickMin, 
% and CtrlVar.ResetThicknessToMinThickness true, then h is first modified accordingly.)
h=s-b;
[b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);

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


if exist('ub','var')==0   % needed if restart file is from the time when ud was simply u
    ub=zeros(MUA.Nnodes,1) ; vb=zeros(MUA.Nnodes,1) ; ud=zeros(MUA.Nnodes,1) ; vd=zeros(MUA.Nnodes,1) ;
    dubdt=zeros(MUA.Nnodes,1) ; dvbdt=zeros(MUA.Nnodes,1) ; duddt=zeros(MUA.Nnodes,1) ; dvddt=zeros(MUA.Nnodes,1) ;
    dubdtm1=zeros(MUA.Nnodes,1) ; dvbdtm1=zeros(MUA.Nnodes,1) ; duddtm1=zeros(MUA.Nnodes,1) ; dvddtm1=zeros(MUA.Nnodes,1) ;
end



if CtrlVar.ResetTime==1 ;
    time=CtrlVar.time;
    CurrentRunStepNumber=0; 
    fprintf(CtrlVar.fidlog,' Time reset to %-g \n',time);
end

if CtrlVar.ResetTimeStep==1 ;
    dt=CtrlVar.dt;
    fprintf(CtrlVar.fidlog,' Time-step reset to %-g \n',dt);
end

fprintf(CtrlVar.fidlog,' Read restart file %s.  Starting restart run at t=%-g with dt=%-g \n',CtrlVar.NameOfRestartFiletoRead,time,dt);




if  time> CtrlVar.TotalTime
    fprintf(CtrlVar.fidlog,' Time at restart (%-g) larger than total run time (%-g) and run  is terminated. \n',time,CtrlVar.TotalTime) ;
    return
end

MeshChanged=0;
if CtrlVar.ReadInitialMesh==1
    fprintf(CtrlVar.fidlog,' On restart loading an initial mesh from %s \n ',CtrlVar.ReadInitialMeshFileName);
    fprintf(CtrlVar.fidlog,' This new mesh will replace the mesh in restart file. \n');
    
    MUAold=MUA;
    clear MUA
    
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
    
    MeshChanged=HasMeshChanged(MUAold,MUA);
    
    if ~MeshChanged
        fprintf(CtrlVar.fidlog,' The new mesh in %s is found to be identical to the mesh in restart file (%s).\n',CtrlVar.ReadInitialMeshFileName,CtrlVar.NameOfRestartFiletoRead);
        fprintf(CtrlVar.fidlog,' No mapping of variables between meshes needed. \n');
    end
    
end
    
for I=1:CtrlVar.RefineMeshOnRestart
    fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
    MUAold=MUA;
    [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    MeshChanged=1;
end

if CtrlVar.TriNodes~=size(MUA.connectivity,2)
    fprintf(CtrlVar.fidlog,' Changing element type at restart step because element type defined in Ua2D_InitialUserInput (%-i) different from one found in restart file (%-i) \n',CtrlVar.TriNodes,MUA.nod);
    MUAold=MUA;
    [MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,CtrlVar.TriNodes);
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    MeshChanged=1;
end

if MeshChanged
    fprintf(CtrlVar.fidlog,' Grid changed, all variables mapped from old to new grid \n ');
    
    
    OutsideValues=0;
    [s,b,h,S,B,rho,AGlen,n,C,m,GF,ub,vb,ud,vd]=...
        MapQuantitiesToNewFEmesh(CtrlVar,MUA,MUAold,h,CtrlVar.time,OutsideValues,...
        ub,vb,ud,vd);
    
    %     OutsideValues=0;
    %    x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
    %     [s,b,ub,vb,ud,vd]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValues,s,b,ub,vb,ud,vd);
    clear MUAold
    %h=s-b;
    
    %     Nnodes=length(x);
    dubdt=zeros(MUA.Nnodes,1); dvbdt=zeros(MUA.Nnodes,1);
    duddt=zeros(MUA.Nnodes,1); dvddt=zeros(MUA.Nnodes,1);
    dhdt=zeros(MUA.Nnodes,1);
    dsdt=zeros(MUA.Nnodes,1);
    dubdtm1=zeros(MUA.Nnodes,1); dvbdtm1=zeros(MUA.Nnodes,1);
    duddtm1=zeros(MUA.Nnodes,1); dvddtm1=zeros(MUA.Nnodes,1);
    dhdtm1=zeros(MUA.Nnodes,1);
    
    
    l=UaLagrangeVariables;
    %MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    
end



%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);

%% In principle these calls to Define.. routines should not be needed
[~,~,S,B,alpha]=GetGeometry(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,'SB');
if any(isnan(S)) ; error(' S returned by DefineGeometry contains NaN') ; end
if any(isnan(B)) ; error(' B returned by DefineGeometry contains NaN') ; end


[rho,rhow,g]=GetDensities(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B);
rho=rho+zeros(length(MUA.coordinates),1);  % make sure that rho is a nodal vector
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);

[C,m]=GetSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
[AGlen,n]=GetAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
[as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);

if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(C))
    error(' C is element-based but does not have same number of elements as there are elements in mesh ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(C))
    error(' C is node-based but does not have same number of elements as there are nodes in mesh ')
end


BCs=GetBoundaryConditions(CtrlVar.Experiment,CtrlVar,MUA,BCs,CtrlVar.time,s,b,h,S,B,ub,vb,ud,vd,GF);

if CtrlVar.doplots==1 && CtrlVar.PlotBCs==1 ;
    
    figure
    PlotBoundaryConditions(CtrlVar,MUA,BCs);
    
end



end