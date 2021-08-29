function [UserVar,RunInfo,LSF,Mask,l,LSFqx,LSFqy]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%
%    df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
%
%

narginchk(7,8)
nargoutchk(7,7)


persistent nCallCounter

if isempty(nCallCounter)
    nCallCounter=0;
end

if ~CtrlVar.DevelopmentVersion
    
    error('LevelSetEquation:Development','LevelSetEquation is in deveopment. Do not use.')
    
end



if ~CtrlVar.LevelSetMethod
    LSF=F0.LSF;
    l=[];
    Mask=[] ;
    return
end

if nargin<8
    l=[];
end


if CtrlVar.CalvingLaw=="-No Ice Shelves-"
    
    % I think this could be simplified, arguably no need to calculate signed distance
    % in this case. Presumably could just define the LSF as distance from flotation, ie
    % h-hf.
    [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F1.LSF,0);
    Mask=CalcMeshMask(CtrlVar,MUA,LSF,0);
    return
end

if  ~isfield(CtrlVar,'LevelSetPhase') ||   isempty(CtrlVar.LevelSetPhase) || CtrlVar.LevelSetPhase==""
    % So the Level Set Phase was not prescribed in the call,
    if  mod(nCallCounter,CtrlVar.LevelSetInitialisationInterval)==0
        CtrlVar.LevelSetPhase="Initialisation and Propagation and FAB" ;
    else
        CtrlVar.LevelSetPhase="Propagation and FAB" ;
    end
    nCallCounter=nCallCounter+1;
end
%% Initialisation phase
if  contains(CtrlVar.LevelSetPhase,"Initialisation")
    
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
    end
    
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
    %%
    
    
    %% Now use the Fixed-point approach. That is, solve the level-set equation using only the non-linear FAB diffusion term
    CtrlVar.LSF.L=0 ;   % The level-set equation only (i.e. without the pertubation term)
    CtrlVar.LSF.P=1 ;   % P is the pertubation term (i.e. the FAB term)
    CtrlVar.LSF.T=0 ;
    CtrlVar.LevelSetTheta=1;  % Here use backward Euler to ensure that the final level set is not affected by the initial guess
    
    % This step does not advance the solution forward in time, just solves
    % the diffusion term
    [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
    F1.LSF=LSF ; F1.LSFqx=LSFqx; F1.LSFqy=LSFqy;
    F0.LSF=LSF ; F0.LSFqx=LSFqx; F0.LSFqy=LSFqy;
    
    
    if ~RunInfo.LevelSet.SolverConverged || CtrlVar.LevelSetTestString=="-pseudo-forward-"

        
        % If fixed-point solution did not converge, do a pseudo-forward time stepping
        %
        %
        CtrlVar.LSF.T=1 ;CtrlVar.LSF.L=0 ;  CtrlVar.LSF.P=1 ;  % Pseudo-forward, using T and P term (no time update and backward Euler)
        CtrlVar.LevelSetTheta=1;
        N=0;
        
        F1.LSF=LSF ; F0.LSF=LSF ;
        Nmax=200;  dtOld=CtrlVar.dt ; dtFactor=2;
        while true
            N=N+1;
            F0.LSF=F1.LSF ;
            
            
            
            [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
            F1.LSF=LSF;
            
            dtBefore=CtrlVar.dt;
            dtNew=CtrlVar.dt*dtFactor ;
            CtrlVar.dt=min(dtNew,1000*dtOld);
            
            dLSFdtMax=max(F1.LSF-F0.LSF)/CtrlVar.dt ;
            
            
            if N>5 && dLSFdtMax < CtrlVar.LevelSetPseudoForwardTolerance
                break
            end
            
            if N>Nmax
                fprintf("Level set solver did not converge despite repeated atempts. \n")
                fprintf("Returning last iterate. Level-set solution might be inaccurate. \n")
                break
            end
            fprintf("time=%f \t dtNew=%g \t dtBefore=%g \t dt=%f \n",CtrlVar.time,dtNew,dtBefore,CtrlVar.dt)
          
        end
        BCs.LSFFixedNode=LSFFixedNodeUnmodified;
        BCs.LSFFixedValue=LSFFixedValueUnmodified;
    end
end

%% Propagation phase, with or without FAB
if contains(CtrlVar.LevelSetPhase,"Propagation")
    
    CtrlVar.LevelSetTheta=0.5;
    
    if contains(CtrlVar.LevelSetPhase,"Propagation and FAB")
        
        CtrlVar.LSF.L=1 ;
        CtrlVar.LSF.P=1 ;
        CtrlVar.LSF.T=1 ;
    else % Propagation only
        CtrlVar.LSF.L=1 ;   % The level-set equation only (i.e. without the pertubation term)
        CtrlVar.LSF.P=0 ;
        CtrlVar.LSF.T=1 ;
        
    end
    
    
    dtOriginal=CtrlVar.dt ;  tEnd=CtrlVar.time+CtrlVar.dt ;
    dtFactor=2; NitDesired=4;  Ntries=0 ;  NtriesMax=20;
    
    
    
    
    while true
        
        
        [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
        
        if RunInfo.LevelSet.SolverConverged
            
            % OK, it converged, advance solution, update time
            
            F1.LSF=LSF;                              % updating solution
            CtrlVar.time=CtrlVar.time+CtrlVar.dt ;   % advancing time
            
            if CtrlVar.time >= tEnd
                % end of time step reached, break out of loop
                break
            end
        
            F0.LSF=LSF ;  % because I will now be doing another sub-forward step
            
            % OK, the step converged but end time of the current time step has not yet been reached (this will happen if previously the
            % time step needed to be reduced due to loss of convergence).
            % Now selecting a new time step based on numer of NR iterations
            dtBefore=CtrlVar.dt;
            dtNew=CtrlVar.dt*(NitDesired/RunInfo.LevelSet.Iterations(RunInfo.LevelSet.iCount));
            CtrlVar.dt=max(min(dtBefore*dtFactor,dtNew),dtBefore/dtFactor);
            
            if (CtrlVar.time+CtrlVar.dt) > tEnd  % don't overstep
                CtrlVar.dt=tEnd-CtrlVar.time;
            end
            
            
            
        else
            
            Ntries=Ntries+1;
            
            if Ntries>NtriesMax
                
                fprintf("Level set solver did not converge despite repeated atempts. \n")
                fprintf("Returning last iterate. Level-set solution might be inaccurate. \n")
                break
                
            elseif Ntries==1
                
                % before reducing time step, first try backward Euler
                CtrlVar.LevelSetTheta=1;
                dtBefore=CtrlVar.dt;
                dtNew=CtrlVar.dt ;
                
                
            else
                % oops, did not converge, so decrease time step and do not advance solution or time,
                % and try again.
                CtrlVar.LevelSetTheta=1;  % Backward Euler
                dtBefore=CtrlVar.dt;
                dtNew=dtBefore/10 ;
                CtrlVar.dt=dtNew;
                fprintf("Level set solver did not converge. Reducing time step and attempting solve again. \n")
                
            end
            
        end
        
        fprintf("time=%f \t tEnd=%f \t dtNew=%g \t dtOld=%g \t dt=%g  \n",CtrlVar.time,tEnd,dtNew,dtBefore,CtrlVar.dt)
        
        
    end
    
    CtrlVar.dt=dtOriginal ;
    fprintf("LSF: time=%f \t tEnd=%f  \t dt=%g  \n",CtrlVar.time,tEnd,CtrlVar.dt)
    
end


% figure ; yyaxis left ; plot(F1.x,LSF,'.c') ; hold on ; plot(F1.x,F1.LSF,'ob') ; yyaxis right ; plot(F1.x,F1.LSF-LSF,'.r')

if ~RunInfo.LevelSet.SolverConverged
    % oops
    warning('LevelSetEquation:NoConvergence','LSF did not converge')
    fprintf('LevelSetEquation:  Solver did not converge.\n')
    fprintf('LevelSetEquation:  Returning last iterate.\n')
end



Mask=CalcMeshMask(CtrlVar,MUA,LSF,0);



if CtrlVar.LevelSetInfoLevel>=10 && CtrlVar.doplots
    
    F1.LSF=LSF ; % here needed for plotting
    [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1);
    
end


end



