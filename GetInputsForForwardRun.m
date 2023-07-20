function  [UserVar,RunInfo,MUA,BCs,F,l]=GetInputsForForwardRun(UserVar,CtrlVar,RunInfo)
          
narginchk(3,3) 
nargoutchk(6,6)

F=UaFields;
EleSizeScalarField=[]; 

F.time=CtrlVar.time ; 

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
   
    [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,F,EleSizeScalarField);
    
    
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

PrintInfoAboutElementsSizes(CtrlVar,MUA); 

if ~isempty(CtrlVar.SaveInitialMeshFileName)
    save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,' MUA was saved in %s .\n',CtrlVar.SaveInitialMeshFileName);
end

F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ; 

BCs=BoundaryConditions;

if CtrlVar.OnlyMeshDomainAndThenStop
    
    % plot mesh, even if PlotMesh not true
    if  CtrlVar.doplots==1
        CtrlVar.PlotMesh=1;
        FindOrCreateFigure("Mesh")
        PlotMuaMesh(CtrlVar,MUA); 
    end
    
    BCs=[]; F=[] ; l=[] ; 
    
    fprintf(CtrlVar.fidlog,' Exiting because CtrlVar.OnlyMeshDomainAndThenStop set to true. \n');
    return
end





[UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,"-s-b-S-B-rho-rhow-g-");

[UserVar,F]=GetCalving(UserVar,CtrlVar,MUA,F,BCs);  % Level Set

if CtrlVar.LevelSetMethod
    % Now this is a special initial case right at the beginning of the run
    % where LSF has been defined ahead of any uv or uvh solutions.
    % Here set all ice thicknesses strickly downstream of the zero level of the LSF to min.
    fprintf("Setting ice thicknesses downstream of calving fronts to the minimum prescribed value of %f .\n",CtrlVar.LevelSetMinIceThickness)
    F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    F.h(F.LSFMask.NodesOut)=CtrlVar.LevelSetMinIceThickness;
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
end


[UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F);
[UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);

F=StartVelocity(CtrlVar,MUA,BCs,F);  % initialize 

[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F);

F=StartVelocity(CtrlVar,MUA,BCs,F);  % modify based on BCs

[UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,BCs,F);

[UserVar,F]=GetStartVelValues(UserVar,CtrlVar,MUA,BCs,F);


l=UaLagrangeVariables; 


F.dubdt=zeros(MUA.Nnodes,1) ;
F.dvbdt=zeros(MUA.Nnodes,1) ;
F.duddt=zeros(MUA.Nnodes,1) ;
F.dvddt=zeros(MUA.Nnodes,1) ;
F.dhdt=zeros(MUA.Nnodes,1) ; 


% is prescribed initial thickness distribution consistent the thickness boundary conditions?
MLC=BCs2MLC(CtrlVar,MUA,BCs);

if ~isempty(MLC.hL)
    BCsThicknessError=norm(MLC.hL*F.h-MLC.hRhs)/norm(MLC.hRhs+1000*eps) ;
    
    if BCsThicknessError > 0.1
        
        fprintf('\n The user-defined initial ice thickness distribution is inconsistent with the user-defined thickness boundary conditions.\n')
        fprintf(' Node \t \t h (defined) \t \t \t h (BCs) \n')
        fprintf('%i : \t \t %f \t \t %f \n ',[BCs.hFixedNode  MLC.hL*F.h MLC.hRhs]')
        
        fprintf('Redefine either the initial ice thickness distribution or the boundary conditions for the ice thickness. \n')
        error('GetInputsForForwardRun:IncorrectUserInputs','User inputs are inconsistent')
        
    end
end

%[UserVar,F]=DefineInputFieldsModifications(UserVar,CtrlVar,F); 


end



