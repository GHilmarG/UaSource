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
    
end
%% Initialisation phase
if  contains(CtrlVar.LevelSetPhase,"Initialisation")
    
    
    F0.LSF=F1.LSF ;  % update F0
    Threshold=0 ;    % Level Set value
    % Here F0.LSF is the original, and F1.LSF will be the re-initilized LSF
    % fix the LSF field for all nodes of elements around the level.
    if CtrlVar.LevelSetInitBCsZeroLevel
        % Use BCs to fix the level set over all elements that the level
        % goes through. This ensures that the level can not shift during
        % initialisation.
        Mask=CalcMeshMask(CtrlVar,MUA,F0.LSF,Threshold);
        BCs.LSFFixedNode=[BCs.LSFFixedNode ; find(Mask.NodesOn)];  % add the nodes of the "On" elements, ie all elements containing the zero level
        BCs.LSFFixedValue=[BCs.LSFFixedValue ; F0.LSF(Mask.NodesOn) ];
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
    
    
    [LSF,UserVar,RunInfo]=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,F0.LSF,xC,yC);
    F0.LSF=LSF ;
    F1.LSF=LSF ;
    %%
    
    
    %% Now use the Fixed-point approach. That is, solve the level-set equation using only the non-linear FAB diffusion term
    CtrlVar.LSF.L=0 ;   % The level-set equation only (i.e. without the pertubation term)
    CtrlVar.LSF.P=1 ;   % P is the pertubation term (i.e. the FAB term)
    CtrlVar.LSF.T=0 ;
    CtrlVar.LevelSetTheta=1;  % Here use backward Euler to ensure that the final level set is not affected by the initial guess
    
    [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
    F1.LSF=LSF ; F1.LSFqx=LSFqx; F1.LSFqy=LSFqy;
    F0.LSF=LSF ; F0.LSFqx=LSFqx; F0.LSFqy=LSFqy;
    
    
    if ~RunInfo.LevelSet.SolverConverged || CtrlVar.LevelSetTestString=="-pseudo-forward-"
        
        % If fixed-point solution did not converge, do a pseudo-forward time stepping
        %
        % But, the fix-point approach works fine and I don't think this
        % is ever needed, keep it in here just in case.
        %
        CtrlVar.LSF.T=1 ;CtrlVar.LSF.L=0 ;  CtrlVar.LSF.P=1 ;
        CtrlVar.LevelSetTheta=1;
        N=0; fprintf("N:%i norm(F1.LSF-F0.LSF)/norm(F0.LSF)=%g \n ",N,norm(F1.LSF-F0.LSF)/norm(F0.LSF))
        F1.LSF=LSF ; F0.LSF=LSF ;
        Nmax=100; tol=1e-4; factor=2;  dtOld=CtrlVar.dt ;
        while true
            N=N+1;
            F0.LSF=F1.LSF ;
            CtrlVar.dt=min([CtrlVar.dt*factor,dtOld*1000]);
            [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
            F1.LSF=LSF;
            
            dlsf=norm(F1.LSF-F0.LSF)/norm(F0.LSF);
            fprintf("N:%i norm(F1.LSF-F0.LSF)/norm(F0.LSF)=%g \n ",N,norm(F1.LSF-F0.LSF)/norm(F0.LSF))
            
            if dlsf < tol || N>Nmax
                break
            end
            
        end
        CtrlVar.dt=dtOld;
        F1.LSF=LSF ; F0.LSF=LSF ;
        %%
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
    
    [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
    
end




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



