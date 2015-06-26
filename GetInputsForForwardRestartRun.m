function    [CtrlVar,MUA,BCs,time,dt,s,b,S,B,ub,vb,ud,vd,l,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
    dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
    GF,GLdescriptors]=GetInputsForForwardRestartRun(CtrlVar)



RestartFile=CtrlVar.NameOfRestartFiletoRead;
Experiment=CtrlVar.Experiment;


dtTemp=CtrlVar.dt;
try
    
    load(RestartFile,'CtrlVarInRestartFile','MUA','BCs','time','dt','s','b','S','B','h','ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
        'Itime','a0','da0dt','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
        'GLdescriptors','l');
    
catch exception
    fprintf(CtrlVar.fidlog,'%s \n',exception.message);
    error('could not load restart file ')
end

if exist('MUA','var')==0
    fprintf(' The variable MUA not found in restart file. Try to read connectivity and coordinates from restart file and then to create MUA \n')
    load(RestartFile,'connectivity','coordinates')
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

CtrlVar.time=time; CtrlVar.dt=dt;

if CtrlVar.ResetTime==1 ;
    time=0;  CtrlVar.time=time;
    Itime=0;
    fprintf(CtrlVar.fidlog,' Time reset to %-g \n',time);
end

if CtrlVar.ResetTimeStep==1 ;
    dt=dtTemp ;  CtrlVar.dt=dt;
    fprintf(CtrlVar.fidlog,' Time-step reset to %-g \n',dt);
end

fprintf(CtrlVar.fidlog,' Read restart file %s.  Starting restart run at t=%-g with dt=%-g \n',RestartFile,time,dt);

CtrlVar.RestartTime=time;


if  time> CtrlVar.TotalTime
    fprintf(CtrlVar.fidlog,' Time at restart (%-g) larger than total run time (%-g) and run  is terminated. \n',time,CtrlVar.TotalTime) ;
    return
end

CtrlVar.MeshChanged=0;
if CtrlVar.ReadInitialMesh==1
    fprintf(CtrlVar.fidlog,' On restart loading an initial mesh from %s \n ',CtrlVar.ReadInitialMeshFileName);
    fprintf(CtrlVar.fidlog,' This new mesh will replace the mesh in restart file. \n')
    CtrlVar.ReadInitialMesh=0;
    MUAold=MUA;
    try
        load(CtrlVar.ReadInitialMeshFileName,'MUA')
        MUA=UpdateMUA(CtrlVar,MUA);
    catch
        load(CtrlVar.ReadInitialMeshFileName,'coordinates','connectivity')
        MUA=CreateMUA(CtrlVar,connectivity,coordinates);
        clear connectivity coordinates
    end
    CtrlVar.MeshChanged=1;
end
    
for I=1:CtrlVar.RefineMeshOnRestart
    fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
    MUAold=MUA;
    [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    CtrlVar.MeshChanged=1;
end

if CtrlVar.TriNodes~=size(MUA.connectivity,2)
    fprintf(CtrlVar.fidlog,' Changing element type at restart step because element type defined in Ua2D_InitialUserInput (%-i) different from one found in restart file (%-i) \n',CtrlVar.TriNodes,MUA.nod);
    MUAold=MUA;
    [MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,CtrlVar.TriNodes);
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    CtrlVar.MeshChanged=1;
end

if CtrlVar.MeshChanged
    fprintf(CtrlVar.fidlog,' Grid changed, all variables mapped from old to new grid \n ');
    
    
    OutsideValues=0;
    [s,b,h,S,B,rho,AGlen,n,C,m,GF,ub,vb,ud,vd]=...
        MapQuantitiesToNewFEmesh(CtrlVar,MUA,MUAold,h,time,OutsideValues,...
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
    
    
    %MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
    
end



%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);

%% In principle these calls to Define.. routines should not be needed
[~,~,S,B,alpha]=DefineGeometry(Experiment,CtrlVar,MUA,time,'SB');
if any(isnan(S)) ; error(' S returned by DefineGeometry contains NaN') ; end
if any(isnan(B)) ; error(' B returned by DefineGeometry contains NaN') ; end


[rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);
rho=rho+zeros(length(MUA.coordinates),1);  % make sure that rho is a nodal vector
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);

[C,m]=GetSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
[AGlen,n]=GetAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
[as,ab]=GetMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);

if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(C))
    error(' C is element-based but does not have same number of elements as there are elements in mesh ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(C))
    error(' C is node-based but does not have same number of elements as there are nodes in mesh ')
end


BCs=GetBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF);

if CtrlVar.doplots==1 && CtrlVar.PlotBCs==1 ;
    
    figure
    PlotBoundaryConditions(CtrlVar,MUA,BCs);
    
end

CtrlVar.MeshChanged=0;

end