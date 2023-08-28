function [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%
% 
% 
% $$\partial_t f +  \mathbf{v} \cdot \nabla f  - \nabla \cdot (\kappa \nabla f) = c \, \|(\nabla f)\|$$
% 
%
%
%    df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
%
%

narginchk(7,8)
nargoutchk(6,6)


persistent nCallCounter

if isempty(nCallCounter)
    nCallCounter=0;
end


if CtrlVar.LevelSetMethod && ~CtrlVar.DevelopmentVersion
    
    error('LevelSetEquation:Development','LevelSetEquation is in deveopment. Do not use.')
    
end




if ~CtrlVar.LevelSetMethod
 
    LSF=[];
    l=[];
    LSFqx=[];
    LSFqy=[]; 
    return
end

if any(isnan(F0.c))
    fprintf("LevelSetEquationSolver: Level set is not evolved because calving rate (c) contains nan. \n")
    LSF=F1.LSF;
    Mask=[];
    l=[];
    LSFqx=[];
    LSFqy=[]; 
    return
end


if nargin<8
    l=[];
end






CtrlVar.LevelSetTheta=0.5;
% [rP,rL,rTP,rTL,rPTL]=CalcLSFconstfunctionTLPterms(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l,F1.LSF);
% if rP/rL> 1e-3
%     fprintf("LevelSetEquationSolver:Automated re-initialisation. \n")
%     
%     CtrlVar.LevelSetFixPointSolverApproach="-PTS-FFP-";
%     CtrlVar.LevelSetReinitializePDist=1;
%     CtrlVar.LevelSetPseudoFixPointSolverMaxIterations=10;
%     CtrlVar.LevelSetTheta=1;
%     [UserVar,RunInfo,LSFini,Mask,l,~,~,BCsLevel]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
%     F0.LSF=LSFini;
%     F1.LSF=LSFini;
% 
% 
%     CalcLSFconstfunctionTLPterms(UserVar,RunInfo,CtrlVar,MUA,BCsLevel,F0,F1,l,LSFini);
% 
%     fprintf("LevelSetEquationSolver:After re-initialisation: \n")
%     % and then recalculate as it would be done
%     CtrlVar.LevelSetTheta=0.5 ;
%     CalcLSFconstfunctionTLPterms(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l,F1.LSF);
% end


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

%CtrlVar.LineUpGLs=false ; Threshold=0 ; 
%[xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Threshold);
%[SignedDist,UserVar,RunInfo]=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,F1.LSF,xc,yc);
% Slope=abs(F1.LSF./SignedDist);
% DistMin=50e3 ; MinSlope=min(Slope(abs(SignedDist)>DistMin));
% fprintf("\n\n =======================  min(Slope)=%f \n\n",MinSlope)

if  contains(CtrlVar.LevelSetPhase,"Initialisation")
  
     [UserVar,RunInfo,LSF,l,LSFqx,LSFqy,BCs]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
     F0.LSF=LSF ; F1.LSF=LSF ;

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
                
                fprintf("LevelSetEquationSolver: Level set solver did not converge despite repeated atempts. \n")
                fprintf("Returning last iterate. Level-set solution might be inaccurate. \n")
                break
                
            elseif Ntries==1
                
                % before reducing time step, first try backward Euler
                CtrlVar.LevelSetTheta=1;
                dtBefore=CtrlVar.dt;
                dtNew=CtrlVar.dt ;
                fprintf("LevelSetEquationSolver: Level set solver did not converge. Trying backward Euler. \n")

            elseif Ntries==2
                
                fprintf("LevelSetEquationSolver: Level set solver did not converge. Performing a new re-initialisation \n")
                CtrlVar.LevelSetReinitializePDist=false ; 
                [UserVar,RunInfo,LSF,Mask,l,LSFqx,LSFqy]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
                F0.LSF=LSF ; F1.LSF=LSF ;
                
            else
                % oops, did not converge, so decrease time step and do not advance solution or time,
                % and try again.
                CtrlVar.LevelSetTheta=1;  % Backward Euler
                dtBefore=CtrlVar.dt;
                dtNew=dtBefore/10 ;
                CtrlVar.dt=dtNew;
                fprintf("LevelSetEquationSolver: Level set solver did not converge. Reducing time step and attempting solve again. \n")
                
            end
            
        end
        
        fprintf("time=%f \t tEnd=%f \t dtNew=%g \t dtOld=%g \t dt=%g  \n",CtrlVar.time,tEnd,dtNew,dtBefore,CtrlVar.dt)
        
        
    end
    
    CtrlVar.dt=dtOriginal ;
    fprintf("LevelSetEquationSolver: LSF time=%f \t tEnd=%f  \t dt=%g  \n",CtrlVar.time,tEnd,CtrlVar.dt)
    
end


% figure ; yyaxis left ; plot(F1.x,LSF,'.c') ; hold on ; plot(F1.x,F1.LSF,'ob') ; yyaxis right ; plot(F1.x,F1.LSF-LSF,'.r')

if ~RunInfo.LevelSet.SolverConverged
    % oops
    warning('LevelSetEquation:NoConvergence','LSF did not converge')
    fprintf('LevelSetEquation:  Solver did not converge.\n')
    fprintf('LevelSetEquation:  Returning last iterate.\n')
end







if CtrlVar.LevelSetInfoLevel>=10 && CtrlVar.doplots
    
    F1.LSF=LSF ; % here needed for plotting
    [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1);
    
end


end



