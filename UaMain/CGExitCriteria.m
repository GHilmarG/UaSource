function [Exit,ExitInfo]=CGExitCriteria(CGExitOptions,ExitInfo,Iteration,J,sGs,cgInfo,lsInfo,Misfit)

%%
%
%   [Exit,ExitInfo]=CGExitCriteria(CtrlVar,ExitInfo,Iteration,J,sGs,cgInfo,lsInfo,Misfit)
%
% Exit criteria for a non-linear conjugated-gradient optimisation using
% NewConjugatedGradMetric.m and LineSearchWolfe.m .
%
% Called once per iteration, after J and the Riesz-mapped gradient have been
% updated. Returns Exit=true when the optimisation should stop, and may instead
% request a CG restart through ExitInfo.ForceRestart .
%
%
%% Inputs
%
%   CtrlVar     options, see below
%   ExitInfo    state, carried between calls. Pass [] on the first call.
%   Iteration   iteration counter, 0 for the initial call before any step
%   J           current value of the cost function
%   sGs         ||s||_G^2 = s'*G*s , with s=-dJdp the Riesz-mapped steepest
%               descent direction. NOT norm(dJdp)^2 , see the note below.
%   cgInfo      as returned by NewConjugatedGradMetric. [] on the first call, or
%               whenever the direction was set to steepest descent by the caller
%               rather than by that routine.
%   lsInfo      as returned by LineSearchWolfe. [] on the first call.
%   Misfit      (optional) data misfit term, for the discrepancy principle. Only
%               needed if CtrlVar.Inverse.TargetMisfit is set.
%
% Outputs
%
%   Exit        true if the optimisation should stop
%   ExitInfo    state, with
%
%                 ExitInfo.Flag          see the list of exit flags below
%                 ExitInfo.Message       text describing Flag
%                 ExitInfo.ForceRestart  true if the caller should discard the
%                                        CG direction and use steepest descent
%                                        on the next iteration
%                 ExitInfo.Decrement     current decrement estimate, 0.5*sGs
%                 ExitInfo.dJ            decrease in J over the last iteration
%                 ExitInfo.nStagnate     consecutive stagnating iterations
%                 ExitInfo.nFuncTotal    running total of Func evaluations
%                 ExitInfo.nGradTotal    running total of gradient evaluations
%                 ExitInfo.History       per-iteration record, see below
%
%
%% Why sGs and not norm(dJdp)
%
% Two reasons, and both matter here.
%
% First, ||s||_G^2 = g' inv(G) g is the squared dual norm of the gradient. The
% l2 norm of the nodal gradient vector is mesh dependent, being a Riemann sum
% with its quadrature weights left out, so a tolerance tuned on one mesh means
% something different on another. The dual norm is mesh independent.
%
% Second, and more useful, it is a decrement. For a quadratic with Hessian H,
%
%       J - J* = 0.5 g' inv(H) g
%
% which is the Newton decrement. Replacing H by the metric G gives
%
%       J - J* ~ 0.5 ||s||_G^2 = ExitInfo.Decrement
%
% so there is an estimate of the remaining decrease in J, in the units of J, and
% without ever forming a Hessian. That is what makes a relative tolerance on it
% meaningful, which is not true of a bare gradient norm.
%
% The estimate is only as good as G is a proxy for H. G is a regularisation-type
% metric, so where the data misfit dominates G is smaller than H and the
% decrement overestimates what is left to gain. Treat it as an order of
% magnitude, not a bound.
%
%
%% Options, all fields of CtrlVar.Inverse, set to defaults if absent or empty
%
%   DecrementTolerance      1e-8    stop when Decrement < DecrementTolerance*|J|
%   DecrementAbsTolerance   0       absolute floor on the same test
%   dJTolerance             1e-10   an iteration counts as stagnating when the
%                                   relative decrease in J falls below this
%   nStagnateMax            3       stop after this many consecutive stagnating
%                                   iterations. The first one triggers a CG
%                                   restart rather than an exit.
%   MaxIterations           100
%   MaxFuncEvaluations      inf     budget, counted across all line searches
%   MaxGradEvaluations      inf
%   TargetMisfit            []      discrepancy principle. If set, stop once
%                                   Misfit <= TargetMisfit.
%   InfoLevel               0
%
%
%% Exit flags
%
%    0   continue
%    1   converged, decrement below tolerance
%    2   stagnated, no useful decrease in J over nStagnateMax iterations
%    3   data misfit reached the target, discrepancy principle
%    4   line search failed on a steepest-descent direction
%    5   maximum number of iterations reached
%    6   Func evaluation budget exhausted
%    7   gradient evaluation budget exhausted
%    8   J or sGs not finite
%   -1   line search reports a non-descent direction
%
% Flags 1 and 3 are convergence. Flags 2 and 4 are the usual way a run on real
% data actually ends. Flag -1 should be unreachable given the sufficient-descent
% safeguard in NewConjugatedGradMetric, and if it appears it points to an
% inconsistency between J and its gradient rather than to anything wrong with
% the optimiser.
%
%
%% How the line-search exit flags are used
%
% A single line-search failure is usually a badly scaled gammaStart rather than
% convergence, so it is not by itself a reason to stop. The rule applied here is
%
%   failure on a CG direction               request a restart and carry on
%   failure on a steepest-descent direction stop, flag 4
%
% the reasoning being that if no progress can be made downhill along -inv(G)*g
% itself, there is no progress left to be had. LineSearchWolfe iExit=4, meaning
% gammaMax was reached, is not treated as a failure, but is counted in
% ExitInfo.nGammaMax since a run where it happens often is one where gammaMax is
% set too small.
%
% iExit=3, a collapsed bracket, usually means J is noisy at the level of
% accuracy being asked for, most often through incompletely converged forward
% solves. It is counted separately in ExitInfo.nNoisy. If it recurs, the fix is
% a tighter forward solve or a looser tolerance here, not a change to the line
% search.
%
%%

if nargin<8 ; Misfit=[] ; end
if nargin<7 ; lsInfo=[] ; end
if nargin<6 ; cgInfo=[] ; end

%% options

CGExitOptions=SetDefault(CGExitOptions,'DecrementTolerance',1e-8) ;
CGExitOptions=SetDefault(CGExitOptions,'DecrementAbsTolerance',0) ;
CGExitOptions=SetDefault(CGExitOptions,'dJTolerance',1e-10) ;
CGExitOptions=SetDefault(CGExitOptions,'nStagnateMax',3) ;
CGExitOptions=SetDefault(CGExitOptions,'MaxIterations',100) ;
CGExitOptions=SetDefault(CGExitOptions,'MaxFuncEvaluations',inf) ;
CGExitOptions=SetDefault(CGExitOptions,'MaxGradEvaluations',inf) ;
CGExitOptions=SetDefault(CGExitOptions,'TargetMisfit',[]) ;
CGExitOptions=SetDefault(CGExitOptions,'InfoLevel',0) ;


%% state

if isempty(ExitInfo) || ~isstruct(ExitInfo)
    ExitInfo=struct ;
    ExitInfo.JLast=inf ;
    ExitInfo.nStagnate=0 ;
    ExitInfo.nFuncTotal=0 ;
    ExitInfo.nGradTotal=0 ;
    ExitInfo.nLineSearchFail=0 ;
    ExitInfo.nGammaMax=0 ;
    ExitInfo.nNoisy=0 ;
    ExitInfo.nRestart=0 ;
    ExitInfo.History=struct('Iteration',[],'J',[],'Decrement',[],'dJ',[],...
        'nFunc',[],'nGrad',[],'iExitLineSearch',[],'teta',[],...
        'SufficientDescentRatio',[],'PowellRatio',[],'CGCorrectionRatio',[],...
        'ddAngle',[],'nCGUpdate',[]) ;
end

ExitInfo.ForceRestart=false ;
ExitInfo.Flag=0 ;
ExitInfo.Message='' ;

Decrement=0.5*sGs ;
ExitInfo.Decrement=Decrement ;

dJ=ExitInfo.JLast-J ;
if ~isfinite(dJ) ; dJ=nan ; end
ExitInfo.dJ=dJ ;

if ~isempty(lsInfo)
    ExitInfo.nFuncTotal=ExitInfo.nFuncTotal+lsInfo.nFuncEvaluations ;
    ExitInfo.nGradTotal=ExitInfo.nGradTotal+lsInfo.nGradEval ;
end

% Was the last search direction a steepest-descent direction rather than a CG
% direction? The counter inside NewConjugatedGradMetric has gone under more than
% one name, so accept either. Note that a missing field is treated as steepest
% descent, which is the safe default but is also silent: if the field name
% changes again, restarts stop being requested and line-search failures exit
% immediately, with nothing to indicate why.
if isempty(cgInfo) || ~isstruct(cgInfo)
    IsSteepestDescent=true ;
elseif isfield(cgInfo,'NumberOfConjGradUpdatesWithoutReset')
    IsSteepestDescent = cgInfo.NumberOfConjGradUpdatesWithoutReset==0 ;
elseif isfield(cgInfo,'ConjGradUpdate')
    IsSteepestDescent = cgInfo.ConjGradUpdate==0 ;
else
    IsSteepestDescent=true ;
    if CGExitOptions.InfoLevel>=1 && Iteration==1
        warning('CGExitCriteria:NoCGCounter',...
            ['cgInfo has no recognised CG update counter, so every direction is ',...
            'being treated as steepest descent and no restarts will be requested.'])
    end
end

Scale=max(abs(J),eps) ;

%% record

ExitInfo.History.Iteration(end+1,1)=Iteration ;
ExitInfo.History.J(end+1,1)=J ;
ExitInfo.History.Decrement(end+1,1)=Decrement ;
ExitInfo.History.dJ(end+1,1)=dJ ;
if isempty(lsInfo)
    ExitInfo.History.nFunc(end+1,1)=0 ;
    ExitInfo.History.nGrad(end+1,1)=0 ;
    ExitInfo.History.iExitLineSearch(end+1,1)=nan ;
else
    ExitInfo.History.nFunc(end+1,1)=lsInfo.nFuncEvaluations ;
    ExitInfo.History.nGrad(end+1,1)=lsInfo.nGradEval ;
    ExitInfo.History.iExitLineSearch(end+1,1)=lsInfo.iExit ;
end
if isempty(cgInfo) || ~isfield(cgInfo,'teta')
    ExitInfo.History.teta(end+1,1)=nan ;
else
    ExitInfo.History.teta(end+1,1)=cgInfo.teta ;
end

% CG degradation diagnostics, see NewConjugatedGradMetric. Recorded here so that
% they can be plotted against the update counter later.
ExitInfo.History.SufficientDescentRatio(end+1,1)=GetField(cgInfo,'SufficientDescentRatio') ;
ExitInfo.History.PowellRatio(end+1,1)=GetField(cgInfo,'PowellRatio') ;
ExitInfo.History.CGCorrectionRatio(end+1,1)=GetField(cgInfo,'CGCorrectionRatio') ;
ExitInfo.History.ddAngle(end+1,1)=GetField(cgInfo,'ddAngle') ;
ExitInfo.History.nCGUpdate(end+1,1)=GetField(cgInfo,'NumberOfConjGradUpdatesWithoutReset') ;

%% the tests, in order of priority

Exit=true ;

if ~isfinite(J) || ~isfinite(sGs)
    
    ExitInfo.Flag=8 ;
    ExitInfo.Message='J or sGs is not finite.' ;
    
elseif ~isempty(lsInfo) && lsInfo.iExit==-1
    
    ExitInfo.Flag=-1 ;
    ExitInfo.Message=['line search reports a non-descent direction. Suspect an ',...
        'inconsistency between J and its gradient.'] ;
    
elseif Decrement < max(CGExitOptions.DecrementTolerance*Scale,CGExitOptions.DecrementAbsTolerance)
    
    ExitInfo.Flag=1 ;
    ExitInfo.Message=sprintf('converged, decrement=%g < %g.',...
        Decrement,max(CGExitOptions.DecrementTolerance*Scale,CGExitOptions.DecrementAbsTolerance)) ;
    
elseif ~isempty(CGExitOptions.TargetMisfit) && ~isempty(Misfit) && Misfit<=CGExitOptions.TargetMisfit
    
    ExitInfo.Flag=3 ;
    ExitInfo.Message=sprintf('data misfit %g has reached the target %g.',...
        Misfit,CGExitOptions.TargetMisfit) ;
    
else
    
    Exit=false ;   % nothing terminal so far, now the recoverable cases
    
    % line-search diagnostics
    if ~isempty(lsInfo)
        
        if lsInfo.iExit==4 ; ExitInfo.nGammaMax=ExitInfo.nGammaMax+1 ; end
        if lsInfo.iExit==3 ; ExitInfo.nNoisy=ExitInfo.nNoisy+1 ; end
        
        if any(lsInfo.iExit==[2 3 -2])
            ExitInfo.nLineSearchFail=ExitInfo.nLineSearchFail+1 ;
            if IsSteepestDescent
                Exit=true ;
                ExitInfo.Flag=4 ;
                ExitInfo.Message=sprintf(...
                    ['line search failed on a steepest-descent direction ',...
                    '(iExit=%i, %s).'],lsInfo.iExit,lsInfo.Message) ;
            else
                ExitInfo.ForceRestart=true ;
                ExitInfo.Message=sprintf(...
                    ['line search failed on a CG direction (iExit=%i), ',...
                    'restarting with steepest descent.'],lsInfo.iExit) ;
            end
        end
        
    end
    
    % stagnation
    if ~Exit && Iteration>0
        
        if isfinite(dJ) && dJ < CGExitOptions.dJTolerance*Scale
            
            ExitInfo.nStagnate=ExitInfo.nStagnate+1 ;
            
            if ExitInfo.nStagnate>=CGExitOptions.nStagnateMax
                Exit=true ;
                ExitInfo.Flag=2 ;
                ExitInfo.Message=sprintf(...
                    ['stagnated, relative decrease in J below %g for %i ',...
                    'consecutive iterations.'],CGExitOptions.dJTolerance,ExitInfo.nStagnate) ;
            elseif ~IsSteepestDescent
                % first sign of stagnation on a CG direction: try a restart
                % before giving up on it
                ExitInfo.ForceRestart=true ;
            end
            
        else
            ExitInfo.nStagnate=0 ;
        end
        
    end
    
    % budgets
    if ~Exit && Iteration>=CGExitOptions.MaxIterations
        Exit=true ; ExitInfo.Flag=5 ;
        ExitInfo.Message=sprintf('maximum number of iterations (%i) reached.',CGExitOptions.MaxIterations) ;
    end
    
    if ~Exit && ExitInfo.nFuncTotal>=CGExitOptions.MaxFuncEvaluations
        Exit=true ; ExitInfo.Flag=6 ;
        ExitInfo.Message=sprintf('Func evaluation budget (%g) exhausted.',CGExitOptions.MaxFuncEvaluations) ;
    end
    
    if ~Exit && ExitInfo.nGradTotal>=CGExitOptions.MaxGradEvaluations
        Exit=true ; ExitInfo.Flag=7 ;
        ExitInfo.Message=sprintf('gradient evaluation budget (%g) exhausted.',CGExitOptions.MaxGradEvaluations) ;
    end
    
end

if ExitInfo.ForceRestart ; ExitInfo.nRestart=ExitInfo.nRestart+1 ; end

ExitInfo.JLast=J ;

%% reporting

if CGExitOptions.InfoLevel>=1
    fprintf(' It %-4i J=%-14.8g dJ=%-12.5g decrement=%-12.5g nFunc=%-5i nGrad=%-5i \n',...
        Iteration,J,dJ,Decrement,ExitInfo.nFuncTotal,ExitInfo.nGradTotal)
end

if CGExitOptions.InfoLevel>=1 && ~isempty(ExitInfo.Message)
    fprintf('   CGExitCriteria: %s \n',ExitInfo.Message)
end

end

%%

function CGExitOptions=SetDefault(CGExitOptions,Field,Value)

if ~isfield(CGExitOptions,Field) || isempty(CGExitOptions.(Field))
    CGExitOptions.(Field)=Value ;
end

end

%%

function v=GetField(S,Field)

% value of S.(Field) if present, NaN otherwise

if isstruct(S) && isfield(S,Field) && ~isempty(S.(Field))
    v=S.(Field) ;
else
    v=nan ;
end

end
