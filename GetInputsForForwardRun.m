function  [MeshChanged,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
    dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1=...
    GetInputsForForwardRun(CtrlVar)

% [CtrlVar,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
%             dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
%             GF]=GetInputsForForwardRun(CtrlVar);


narginchk(1,1)
nargoutchk(34,34)

%for I=1:CtrlVar.nInitialRemeshSteps+1

if CtrlVar.ReadInitialMesh==1
    
    fprintf(CtrlVar.fidlog,' loading an initial mesh from %s \n ',CtrlVar.ReadInitialMeshFileName);
      
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
    
    for I=1:CtrlVar.RefineMeshOnStart
        fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
        [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
        MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
       
    end
    
else
   
    MUA=genmesh2d(CtrlVar,CtrlVar.MeshBoundaryCoordinates);
    
    
    iIt=0;
    while MUA.Nele>1.2*CtrlVar.MaxNumberOfElements && iIt<=2;
        iIt=iIt+1;
        % Note: these changes in MeshSize are not returned
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
    MeshChanged=[];
    
    fprintf(CtrlVar.fidlog,' Exiting beacause CtrlVar.OnlyMeshDomainAndThenStop set to true. \n');
    return
end

%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);

[s,b,S,B,alpha]=GetGeometry(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,'sbSB');
TestVariablesReturnedByDefineGeometryForErrors(MUA,s,b,S,B);


h=s-b;
[rho,rhow,g]=GetDensities(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B);
GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);

[C,m]=GetSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
[AGlen,n]=GetAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
[ub,vb,ud,vd]=StartVelocity(CtrlVar,MUA);
[ub,vb,ud,vd]=GetStartVelValues(CtrlVar.Experiment,CtrlVar,MUA,ub,vb,ud,vd,CtrlVar.time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
[as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);



BCs=GetBoundaryConditions(CtrlVar.Experiment,CtrlVar,MUA,BCs,CtrlVar.time,s,b,h,S,B,ub,vb,ud,vd,GF);

if CtrlVar.doplots==1
    if CtrlVar.PlotBCs;
        figure ; PlotBoundaryConditions(CtrlVar,MUA,BCs);
    elseif CtrlVar.PlotMesh
        figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
    end
end


MeshChanged=0;

end



