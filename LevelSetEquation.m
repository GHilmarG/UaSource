function [UserVar,RunInfo,LSF,Mask,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1)
    %%
    %
    %
    %    df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    %
    
    
    narginchk(7,7)
    nargoutchk(4,5)
    
    
    persistent LastResetTime dLSF
    
    if ~CtrlVar.DevelopmentVersion
       
        error('LevelSetEquation:Development','LevelSetEquation is in deveopment. Do not use.')
        
    end
    
    
    
    if ~CtrlVar.LevelSetMethod
        LSF=F0.LSF;
        lambda=[];
        Mask=[] ; 
        return
    end
    
    if CtrlVar.CalvingLaw=="-No Ice Shelves-" 
         
         % I think this could be simplified, arguably no need to calculate signed distance
         % in this case. Presumably could just define the LSF as distance from flotation, ie
         % h-hf. 
        [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F1.LSF,0);
        Mask=CalcMeshMask(CtrlVar,MUA,LSF,0);
        return 
    end
    
    if isempty(LastResetTime)
        LastResetTime=0 ;
    end
    
    
    switch CtrlVar.LevelSetPhase
        case "Initialisation"
            CtrlVar.LSF.L=0 ;   % The level-set equation only (i.e. without the pertubation term)
            CtrlVar.LSF.P=1 ;   % % P is the pertubation term
            CtrlVar.LevelSetTheta=1;  
            CtrlVar.LevelSetEpsilon=0 ;
            
            % First do a rough re-initialisation using signed distance function.
            % After this I then do a full non-linear FAB solve with the level-set fixed 
            % as boundary conditions on the LSF.
            [F1.LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F1.LSF); 
            F0.LSF=F1.LSF ;
            Mask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
            BCs.LSFFixedNode=[BCs.LSFFixedNode ; find(Mask.NodesOn)];   % fix the LSF field for all nodes of elements around the level.
            BCs.LSFFixedValue=[BCs.LSFFixedValue ; F1.LSF(Mask.NodesOn) ];
            
        case "Propagation"
            CtrlVar.LSF.L=1 ;   % The level-set equation only (i.e. without the pertubation term)
            CtrlVar.LSF.P=0 ;
        case "Propagation and FAB"
            CtrlVar.LevelSetTheta=0.5;
            CtrlVar.LSF.L=1 ;
            CtrlVar.LSF.P=1 ;
        otherwise
            error('safd')
    end
    
    
    
    

    
    % This will actually also do Picard unless CtrlVar.LevelSetSolutionMethod="Newton-Raphson" ;
    [UserVar,RunInfo,LSF,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1);
    
    if ~RunInfo.LevelSet.SolverConverged
        % oops
        error('LevelSetEquation:NoConvergence','LSF did not converge')
        fprintf('LevelSetEquation:  Solver did not converge.\n')
        dtKeep=CtrlVar.dt;
        CtrlVar.dt=CtrlVar.dt/10 ;
        F1.LSF=F0.LSF;
        LSF=F0.LSF;
        for iTheta=1:10
            F1.LSF=LSF ;
            [UserVar,RunInfo,LSF,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1);
        end
        CtrlVar.dt=dtKeep;
    end
    
    
    
    Mask=CalcMeshMask(CtrlVar,MUA,LSF,0);
    
    dLSF=LSF-F0.LSF;
    
    if CtrlVar.LevelSetInfoLevel>=100 && CtrlVar.doplots
        
        F1.LSF=LSF ; % here needed for plotting
        [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1);
        
    end
    
    
end



