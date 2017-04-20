function  [UserVar,MUA,BCs,F,l,GF]=GetInputsForForwardRun(UserVar,CtrlVar)



F=UaFields;

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
        error('Ua:GetInputsForForwardRun:ReadInitialMeshFileName','Input file does not contain expected variables')
    end
    clearvars Temp
    
    for I=1:CtrlVar.RefineMeshOnStart
        fprintf(CtrlVar.fidlog,' All triangle elements are subdivided into four triangles \n');
        [MUA.coordinates,MUA.connectivity]=FE2dRefineMesh(MUA.coordinates,MUA.connectivity);
        MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
       
    end
    
else
   
    [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,CtrlVar.MeshBoundaryCoordinates);
    
    
    iIt=0;
    while MUA.Nele>1.2*CtrlVar.MaxNumberOfElements && iIt<=2
        iIt=iIt+1;
        % Note: these changes in MeshSize are not returned
        if numel(CtrlVar.MeshSize)==1
            CtrlVar.MeshSize=1.2*CtrlVar.MeshSize*sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
        else
            CtrlVar.MeshSize(:,3)=1.2*CtrlVar.MeshSize(:,3)*sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
        end
        fprintf(CtrlVar.fidlog,'Nele=%-i > 1.2*NEleMax=%-i , hence desired meshsize is scaled up \n',MUA.Nele,1.2*CtrlVar.MaxNumberOfElements);
        [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,CtrlVar.MeshBoundaryCoordinates);
        fprintf(CtrlVar.fidlog,'new Nele after scale down is %-i  \n',MUA.Nele);
    end
    
    
    
end

PrintInfoAboutElementsSizes(CtrlVar,MUA)

if ~isempty(CtrlVar.SaveInitialMeshFileName)
    save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,' MUA was saved in %s .\n',CtrlVar.SaveInitialMeshFileName);
end



BCs=BoundaryConditions;

if CtrlVar.OnlyMeshDomainAndThenStop
    
    % plot mesh, even if PlotMesh not true
    if  CtrlVar.doplots==1
        CtrlVar.PlotMesh=1;
        figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
    end
    
    BCs=[]; F=[] ; l=[] ; GF=[] ;
    
    fprintf(CtrlVar.fidlog,' Exiting beacause CtrlVar.OnlyMeshDomainAndThenStop set to true. \n');
    return
end

%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);


[UserVar,F.s,F.b,F.S,F.B,F.alpha]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'sbSB');
TestVariablesReturnedByDefineGeometryForErrors(MUA,F.s,F.b,F.S,F.B);

F.h=F.s-F.b;

[UserVar,F]=GetDensities(UserVar,CtrlVar,MUA,F);

[F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);



[UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F,GF);
[UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F,GF);
[UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);

F=StartVelocity(CtrlVar,MUA,BCs,F);  % initialize 

[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F,GF);

F=StartVelocity(CtrlVar,MUA,BCs,F);  % modify based on BCs

[UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,BCs,F,GF);

[UserVar,F]=GetStartVelValues(UserVar,CtrlVar,MUA,BCs,F,GF);


l=UaLagrangeVariables; 


F.dubdt=zeros(MUA.Nnodes,1) ;
F.dvbdt=zeros(MUA.Nnodes,1) ;
F.duddt=zeros(MUA.Nnodes,1) ;
F.dvddt=zeros(MUA.Nnodes,1) ;
F.dhdt=zeros(MUA.Nnodes,1) ; 


%[UserVar,F]=DefineInputFieldsModifications(UserVar,CtrlVar,F); 


end



