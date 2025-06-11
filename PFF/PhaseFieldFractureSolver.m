




function [MUA,BCs,BCsphi,F]=PhaseFieldFractureSolver(UserVar,RunInfo,CtrlVar,MUA,F,BCs)


%%
%
% Solves a phase-field fracture problem. 
%
% Under development. Do not use!  
%
%  Outer solve:   
%
%       Inner solve:
%
%           For Psi0=Psi, refine mesh by repeatedly solving: 
%                  
%           uv=uv(phi)  
%           phi=phi(Psi0)
%
%       end 
%
%       update Psi:
%       Psi=Psi(uv)
%
% end 
%
%%



CtrlVar.Parallel.BuildWorkers=true;
MUA=UpdateMUA(CtrlVar,MUA) ; 

% currently only spatially constant intact A (ie F.AGlen0), and rhoi allowed

isAconstant=all(F.AGlen==F.AGlen(1));
isrhoConstant=all(F.rho==F.rho(1));

if isAconstant && isrhoConstant
    F.AGlen=zeros(MUA.Nnodes,1)+mean(F.AGlen);
    F.rho=zeros(MUA.Nnodes,1) + 920 ;
else

    error("PhaserFieldFractureSolver:AnotConstant","Currently undamaged A and rho must be constant in the phase field solver.")

end

F.AGlen0=F.AGlen(1);
F.rho0=F.rho(1); 
F.b0=F.b ;
F.h0=F.h ;

Gc=CtrlVar.PhaseFieldFracture.Gc ;
l=CtrlVar.PhaseFieldFracture.l ;


if isempty(F.Psi) 
    F.Psi=zeros(MUA.Nnodes,1) ;
end

if isempty(F.phi) 
    F.phi=zeros(MUA.Nnodes,1) ;  % phi=0, undamaged,
end

% phi=1, fully damaged

% Make initial phi feasible
BCsphi=BoundaryConditions;
CtrlVar.BCs="-phi-" ;
fprintf("Getting BCs for the phase field (phi) \n ")
[UserVar,BCsphi]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCsphi,F) ;

F.phi(BCsphi.hFixedNode)=BCsphi.hFixedValue;


nMeshRefinements=CtrlVar.PhaseFieldFracture.MaxMeshRefinements;
nphiUpdates=CtrlVar.PhaseFieldFracture.MaxUpdates ;

iphiUpdate=0; 
Dphi=nan(100,1); DphiCount=0;
EleSizeMin=CtrlVar.MeshSizeMin;

 

CtrlVar.PhaseFieldFracture.iphiUpdate=0;
CtrlVar.PhaseFieldFracture.MeshRefinement=0;

PlotTitle="initial configuration" ;
lm=UaLagrangeVariables ;


%% Initial solve for velocities using the initial phase field 
%  This would typically be at the start of the run where phi=0 everywhere, i.e. undamaged material.
%
F=PPFphi2F(CtrlVar,MUA,F) ;           % Map from phi over to material parameters
[UserVar,RunInfo,F,lm]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,lm) ;   % Now solve for velocities
[F.Psi,e,eInt]=StrainRateEnergy(CtrlVar,MUA,F,F.AGlen0) ;            % Update Psi
PFFPlots(UserVar,CtrlVar,MUA,F,BCs,BCsphi,F.phi,F.Psi,e,PlotTitle) ;

while true % phi "evolution" loop, i.e. here the driving term Psi is updated

    iMeshRefinements=0 ;
    iphiUpdate= iphiUpdate + 1;

    


    while true   % mesh refinement loop, 
        %

        PlotTitle=sprintf("phi loop %i, mesh refinement %i",iphiUpdate,iMeshRefinements) ;

        CtrlVar.BCs="-uv-" ; 
        [UserVar,BCs]= GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F); 
     

        
        %% Solve uv problem: this depends only on phi from the (previous) phase-field solution
        %
        %
   
        F=PPFphi2F(CtrlVar,MUA,F) ;
        [UserVar,RunInfo,F,lm]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,lm) ;
     
        %%
        
        %% Solve phase-field problem: This depends only on the previous uv solution, and not on previous phi or Psi
        %  phi field

        BCsphi=BoundaryConditions;
        CtrlVar.BCs="-phi-" ; 
        [UserVar,BCsphi]= GetBoundaryConditions(UserVar,CtrlVar,MUA,BCsphi,F); 

        [F.Psi,e,eInt]=StrainRateEnergy(CtrlVar,MUA,F,F.AGlen0) ; % Update Psi

        phiLast=F.phi;

        % Now solve the phase-field equation to arrive at a new phi field
        [UserVar,F.phi]=PFFequation(UserVar,CtrlVar,MUA,BCsphi,Gc,l,F.Psi);
    
        F.phi(F.GF.node > 0.5 )=0 ; % Here manually resetting damage to zero over grounded areas

        if ~isfield(CtrlVar.PhaseFieldFracture,"UpdateRatio")
            CtrlVar.PhaseFieldFracture.UpdateRatio=1;
        end
        
        
        dphi=F.phi-phiLast;
        dphiNorm=(dphi'*MUA.M* dphi)/MUA.Area; 
        fprintf("|dphi|=%g \n",dphiNorm)
        DphiCount=DphiCount+1; 
        Dphi(DphiCount)=dphiNorm;

        F.phi= CtrlVar.PhaseFieldFracture.UpdateRatio*F.phi+(1- CtrlVar.PhaseFieldFracture.UpdateRatio)*phiLast;
        
        PFFPlots(UserVar,CtrlVar,MUA,F,BCs,BCsphi,F.phi,F.Psi,e,PlotTitle) ;


        iMeshRefinements=iMeshRefinements+1;

        if iMeshRefinements > nMeshRefinements
            fprintf("Exiting element refinement loop. Max number of refinement iterations reached. \n")
            break
        end

        MUAold=MUA;
        phiEle=Nodes2EleMean(MUAold.connectivity,F.phi) ;
        
        Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity);
        EleSize=sqrt(Tarea);

        % For the time being, the mesh-refinement criterion is hard-wired in the code.
        % 
        ElementsToBeRefined   = ( phiEle > 0.5 )  &  ( EleSize > EleSizeMin ) ;
        ElementsToBeCoarsened = ( phiEle < 0.5 )  ; 


        nEleRefine=numel(find(ElementsToBeRefined)) ; 


        fprintf("number of elements refined %i \n",nEleRefine)

        if nEleRefine==0

            fprintf("Exiting element refinement loop. No elements to be refined. \n")
            break

        end

        
        CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection' ; CtrlVar.InfoLevelAdaptiveMeshing=1;
        [MUAnew,RunInfo]=LocalMeshRefinement(CtrlVar,RunInfo,MUAold,ElementsToBeRefined,ElementsToBeCoarsened) ;
    

        isDefineF=CtrlVar.PhaseFieldFracture.isDefineF ; % I currently use this with the PFF drivers
        if isDefineF

            Fnew=DefineF(UserVar,CtrlVar,MUAnew) ;
            OutsideValues=[] ;
            [RunInfo,Fnew.ub,Fnew.vb,Fnew.phi,Fnew.C,Fnew.B]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,F.ub,F.vb,F.phi,F.C,F.B);

        else

            Fnew=PFFmapF(UserVar,RunInfo,CtrlVar,MUA,MUAnew,F) ;

        end
       
         F=Fnew;
         MUA=MUAnew;
         MUA=UpdateMUA(CtrlVar,MUA) ; 
         lm=UaLagrangeVariables ; % can I update this instead of redefining?

    end


   

    CtrlVar.PhaseFieldFracture.iphiUpdate=iphiUpdate ;  
    

    if iphiUpdate > nphiUpdates
        break
    end

    if dphiNorm < 1e-4  % for the time being this is hardwired
         break
    end

   PlotTitle=sprintf("phi loop %i, mesh refinement %i",iphiUpdate,iMeshRefinements) ;
   [F.Psi,e,eInt]=StrainRateEnergy(CtrlVar,MUA,F,F.AGlen0) ; % Update Psi
   PFFPlots(UserVar,CtrlVar,MUA,F,BCs,BCsphi,F.phi,F.Psi,e,PlotTitle) ;

    drawnow

end



FindOrCreateFigure("Dphi"); 
semilogy(Dphi,'o-r') ; 
xlabel("iterations")
ylabel("$\|\Delta \phi\|$",interpreter="latex")
title("Change in phase field with iteration")


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







































