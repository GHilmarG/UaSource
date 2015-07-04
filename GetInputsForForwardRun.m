function  [CtrlVar,MUA,BCs,time,dt,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
    dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
    GF]=GetInputsForForwardRun(CtrlVar)

narginchk(1,1)
nargoutchk(36,36)


Experiment=CtrlVar.Experiment;
time=CtrlVar.time;
dt=CtrlVar.dt;
%for I=1:CtrlVar.nInitialRemeshSteps+1

if CtrlVar.ReadInitialMesh==1
    fprintf(CtrlVar.fidlog,' loading an initial mesh from %s \n ',CtrlVar.ReadInitialMeshFileName);
    CtrlVar.ReadInitialMesh=0;
    
    try
        load(CtrlVar.ReadInitialMeshFileName,'MUA')
        MUA=UpdateMUA(CtrlVar,MUA);
    catch
        load(CtrlVar.ReadInitialMeshFileName,'coordinates','connectivity')
        MUA=CreateMUA(CtrlVar,connectivity,coordinates);
        clear connectivity coordinates
    end
    
for I=1:CtrlVar.RefineMeshOnStart
        fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
        [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
        MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
        CtrlVar.MeshChanged=1;
    end
    
else
    CtrlVar.MeshChanged=1;
    MUA=genmesh2d(CtrlVar,CtrlVar.MeshBoundaryCoordinates);
    
    
    iIt=0;
    while MUA.Nele>1.2*CtrlVar.MaxNumberOfElements && iIt<=2;
        iIt=iIt+1;
        if numel(CtrlVar.MeshSize)==1
            CtrlVar.MeshSize=1.2*CtrlVar.MeshSize*sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
        else
            CtrlVar.MeshSize(:,3)=1.2*CtrlVar.MeshSize(:,3)*sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
        end
        fprintf(CtrlVar.fidlog,'Nele=%-i > 1.2*NEleMax=%-i , hence desired meshsize is scaled up \n',MUA.Nele,1.2*CtrlVar.MaxNumberOfElements);
        MUA=genmesh2d(CtrlVar,CtrlVar.MeshBoundaryCoordinates);
        fprintf(CtrlVar.fidlog,'new Nele after scale down is %-i  \n',MUA.Nele);
    end
    
    
    
end

PrintInfoAboutElementsSizes(CtrlVar,MUA)

if ~isempty(CtrlVar.SaveInitialMeshFileName)
    save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,' MUA was saved in %s .\n',CtrlVar.SaveInitialMeshFileName);
end

dubdt=zeros(MUA.Nnodes,1) ; dubdtm1=zeros(MUA.Nnodes,1);
dvbdt=zeros(MUA.Nnodes,1) ; dvbdtm1=zeros(MUA.Nnodes,1);
duddt=zeros(MUA.Nnodes,1) ; duddtm1=zeros(MUA.Nnodes,1);
dvddt=zeros(MUA.Nnodes,1) ; dvddtm1=zeros(MUA.Nnodes,1);
dhdt=zeros(MUA.Nnodes,1) ; dhdtm1=zeros(MUA.Nnodes,1);
dsdt=zeros(MUA.Nnodes,1) ; dbdt=zeros(MUA.Nnodes,1);


BCs=BoundaryConditions;

if CtrlVar.OnlyMeshDomainAndThenStop
    
    % plot mesh, even if PlotMesh not true
    if  CtrlVar.doplots==1
        CtrlVar.PlotMesh=1;
        figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
    end
    
    s=[]; b=[] ; S=[] ; B=[] ; rho=[] ; rhow=[]; alpha=[] ; g=[] ;
    C=[] ; m=[] ; AGlen=[] ; n=[] ; ud=[] ; vd=[] ; ub=[] ;vb=[] ;as=[] ; ab=[] ; GF=[];
    
    
    fprintf(CtrlVar.fidlog,' Exiting\n');
    return
end

%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);

[s,b,S,B,alpha]=GetGeometry(Experiment,CtrlVar,MUA,time,'sbSB');
TestVariablesReturnedByDefineGeometryForErrors(MUA,s,b,S,B);


h=s-b;
[rho,rhow,g]=GetDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);
GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);

[C,m]=GetSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
[AGlen,n]=GetAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
[ub,vb,ud,vd]=StartVelocity(CtrlVar,MUA);
[ub,vb,ud,vd]=GetStartVelValues(Experiment,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
[as,ab]=GetMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);


% [ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB,FixedNormalVelocityNode,FixedNormalVelocityValue]=...
%     DefineBCs(Experiment,CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,GF);

BCs=GetBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF);

if CtrlVar.doplots==1
    if CtrlVar.PlotBCs;
        figure ; PlotBoundaryConditions(CtrlVar,MUA,BCs);
    elseif CtrlVar.PlotMesh
        figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
    end
end


%
% if ~MLC.ubvbTies && ~MLC.hTies ;
%     if CtrlVar.InfoLevel>0
%         fprintf(CtrlVar.fidlog,' No uv ties found \n ');
%     end
%     CtrlVar.SymmSolver='EliminateBCsSolveSystemDirectly';
%     CtrlVar.AsymmSolver='EliminateBCsSolveSystemDirectly';
% end
% %]

%% test if variables are OK
if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(C))
    save TestSave ;
    error(' C is element-based but does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(C))
    save TestSave ;
    error(' C is node-based but does not have same number of elements as there are nodes in mesh. All variables saved in TestSave.mat ')
    
end

if CtrlVar.AGlenisElementBased  && ~(length(MUA.connectivity)==length(AGlen))
    save TestSave ;
    error(' AGlen is element-based but does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.AGlenisElementBased && ~(length(MUA.coordinates) == length(AGlen))
    save TestSave ;
    error(' AGlen is node-based but does not have same number of elements as there are nodes in mesh. All variables saved in TestSAve.mat ')
end

CtrlVar.MeshChanged=0;

end



