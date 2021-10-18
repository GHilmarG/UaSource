function [UserVar,RunInfo,LSF,Mask,l,LSFqx,LSFqy]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)


LSFqx=[] ; LSFqy=[] ; 

% Don't redefine F0.LSF as F1.LSF, doing so would push the solution back in tiem
F1.LSF=F0.LSF ;
Threshold=0 ;    % Level Set value
% Here F0.LSF is the original, and F1.LSF will be the re-initilized LSF
% fix the LSF field for all nodes of elements around the level.
if CtrlVar.LevelSetInitBCsZeroLevel
    % Use BCs to fix the level set over all elements that the level
    % goes through. This ensures that the level can not shift during
    % initialisation.
    Mask=CalcMeshMask(CtrlVar,MUA,F0.LSF,Threshold);
    
    LSFFixedNodeUnmodified=BCs.LSFFixedNode ;
    LSFFixedValueUnmodified=BCs.LSFFixedValue ;
    
    BCs.LSFFixedNode= [LSFFixedNodeUnmodified ; find(Mask.NodesOn)];  % add the nodes of the "On" elements, ie all elements containing the zero level
    BCs.LSFFixedValue=[LSFFixedValueUnmodified ; F0.LSF(Mask.NodesOn) ];
   % figure ;  lgd=PlotBoundaryConditions(CtrlVar,MUA,BCs) ;
end


if  CtrlVar.LevelSetReinitializePDist

    %% After having located the 0 level, now do a rough re-initialisation using signed distance function. After this I then do a full
    % non-linear FAB solve with the level-set fixed as boundary conditions on the LSF.
    % This will in most cases not be needed, but

    if  isfield(CtrlVar,'CtrlVar.LevelSetTestString') &&  contains(CtrlVar.LevelSetTestString,"-xc/yc nodes-")
        xC=F0.x(Mask.NodesOn ) ; yC=F0.y(Mask.NodesOn) ;
    else
        CtrlVar.LineUpGLs=false ;
        [xC,yC]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Threshold);
    end

    % It should be OK to do this with LSF at both 0 and 1 as I have already
    % found the location of the level set for F0 and this will be enforced
    % throught the BCs.
    [LSF,UserVar,RunInfo]=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,F0.LSF,xC,yC);
    F0.LSF=LSF ;
    F1.LSF=LSF ;
end
%%


[UserVar,RunInfo,LSF,l]=LevelSetFixPointSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l) ; 
% var her
return


%% Now use the Fixed-point approach. That is, solve the level-set equation using only the non-linear FAB diffusion term
CtrlVar.LSF.L=0 ;   % The level-set equation only (i.e. without the pertubation term)
CtrlVar.LSF.P=1 ;   % P is the pertubation term (i.e. the FAB term)
CtrlVar.LSF.T=0 ;
CtrlVar.LevelSetTheta=1;  % Here use backward Euler to ensure that the final level set is not affected by the initial guess

% This step does not advance the solution forward in time, just solves
% the diffusion term
l=[] ;  % Here I'm solving a different system potentially, so don't use previous l values
[UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
F1.LSF=LSF ; % F1.LSFqx=LSFqx; F1.LSFqy=LSFqy;
F0.LSF=LSF ; % F0.LSFqx=LSFqx; F0.LSFqy=LSFqy;


if ~RunInfo.LevelSet.SolverConverged || CtrlVar.LevelSetTestString=="-pseudo-forward-"
    
    
    % If fixed-point solution did not converge, do a pseudo-forward time stepping
    %
    %
    CtrlVar.LSF.T=1 ;CtrlVar.LSF.L=0 ;  CtrlVar.LSF.P=1 ;  % Pseudo-forward, using T and P term (no time update and backward Euler)
    CtrlVar.LevelSetTheta=1;
    N=0;
    
    F1.LSF=LSF ; F0.LSF=LSF ;
    Nmax=200;  dtOld=CtrlVar.dt ; dtFactor=1.5;  dtMax=100*dtOld ; l=[] ; 
    while true
        N=N+1;
        F0.LSF=F1.LSF ;
        
        
        
        [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
        F1.LSF=LSF;
        
        dtBefore=CtrlVar.dt;
        dtNew=CtrlVar.dt*dtFactor ;
        CtrlVar.dt=min(dtNew,dtMax);
        
        dLSFdtMax=max(F1.LSF-F0.LSF)/CtrlVar.dt ;
        
        
        if N>5 && dLSFdtMax < CtrlVar.LevelSetPseudoForwardTolerance
            break
        end

        if N > 20  && ~RunInfo.LevelSet.SolverConverged
            CtrlVar.LevelSetTheta=1;  % I would like to use backward Euler, but this occasionally does not converge

        end

        if N>Nmax
            fprintf("Level set solver did not converge despite repeated atempts. \n")
            fprintf("Returning last iterate. Level-set solution might be inaccurate. \n")
            break
        end
        fprintf("time=%f \t dtNew=%g \t dtBefore=%g \t dt=%f \t max rate-of-change=%g \n",CtrlVar.time,dtNew,dtBefore,CtrlVar.dt,dLSFdtMax)
        
    end
    BCs.LSFFixedNode=LSFFixedNodeUnmodified;
    BCs.LSFFixedValue=LSFFixedValueUnmodified;
end

end