function  [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%%
% This is basically a non-linear conjugated-gradient minimizer.
%
% It does a reasonably good job. Importantly it does allow for an arbitrary metric, which here is defined by the metric
% matrix G.  This is something that most optimization packages appear not to allow for.
%
% The line search is done using LineSearchWolfe, which returns a minimum satisfying both Wolfe conditions, i.e. both the
% Armijo rule and the curvature condition.
%
% func is the function to be minimized, 
% 
%   [J,dJdp]=func(p+gamma*d) 
% 
%  where p is the parameter set, gamma is the (scalar) line search step size, and d the search direction. 
%
%
%% Cost-function evaluations
%
% Every evaluation of func within an iteration goes through the nested function PhiCached, which stores gamma, J, dJdp
% and the JGH outputs for each trial step. Two evaluations per iteration are saved as a result:
%
%   1) the accepted point has already been evaluated inside the line search, so J, dJdp and fOuts at the new p are read
%      back from the cache instead of calling func(p) again. That call was a full forward AND adjoint solve.
%
%   2) the slope at gammaStart is passed to LineSearchWolfe through info.slopeStart, so that when gammaStart already
%      satisfies both Wolfe conditions the line search returns without evaluating func at all.
%
% The cache is only valid while p and d are unchanged, so it is cleared at the end of each iteration. It is also
% bypassed if kk_proj actually alters p, since the cached values then refer to a different point.
%
% Note that func is called with four output arguments. This does NOT trigger a Hessian assembly: the anonymous function
% is built in InvertForModelParameters with CtrlVar.JGH.CalcHessian=false for all non-Hessian methods, and that value is
% frozen into the closure.
%
%%

narginchk(8,8)
nargoutchk(3,3)


%% Set various options as based on CtrlVar

if isempty(CtrlVar)

    CtrlVar.Inverse.UaConjugatedGradients.Armijo=1e-4;
    CtrlVar.Inverse.UaConjugatedGradients.WolfeCurvature=0.1;
    CtrlVar.Inverse.UaConjugatedGradients.MaxFuncEvaluationsInLineSearch=15;
    CtrlVar.Inverse.UaConjugatedGradients.InfoLevel=0;

    CtrlVar.Inverse.UaConjugatedGradients.UpdateMethod="-ConjGrad-" ; %{'SteepestDecent','ConjGrad'}
   
    CtrlVar.Inverse.UaConjugatedGradients.Update="HS";
    CtrlVar.Inverse.UaConjugatedGradients.SufficientDescent=0.1;
 
  
    CtrlVar.Inverse.UaConjugatedGradients.DecrementAbsTolerance=1e-10;
    CtrlVar.Inverse.UaConjugatedGradients.dJTolerance=1e-10;
end

% Options for line search, used below by LineSearchWolfe.m
LineSearchOptions.c1=CtrlVar.Inverse.UaConjugatedGradients.Armijo;
LineSearchOptions.c2=CtrlVar.Inverse.UaConjugatedGradients.WolfeCurvature;
LineSearchOptions.MaxFuncEvaluations=CtrlVar.Inverse.UaConjugatedGradients.MaxFuncEvaluationsInLineSearch;
LineSearchOptions.InfoLevel=CtrlVar.Inverse.UaConjugatedGradients.InfoLevel;



if isfield(CtrlVar,"ConjugatedGradientsUpdate")

    fprintf("'CtrlVar.ConjugatedGradientsUpdate' no longer used. \n")
    fprintf("Use insted: 'CtrlVar.Inverse.UaConjugatedGradients.Update' \n")
    error("FieldNoLongerUsed")

end

if isfield(CtrlVar,"ConjugatedGradientsSufficientDescent")
    error("FieldNoLongerUsed")
end

cgInfo.NumberOfConjGradUpdatesWithoutReset=0;
cgInfo.ConjGradAngle=nan;
cgInfo.ddAngle=nan;
cgInfo.teta=nan;
cgInfo.SufficientDescentRatio=nan;
cgInfo.PowellRatio=nan;
cgInfo.CGCorrectionRatio=nan;

% options related to exit criteria, used by CGExitCriteria.m
CGExitOptions.MaxIterations=CtrlVar.Inverse.Iterations;
CGExitOptions.DecrementTolerance=CtrlVar.Inverse.UaConjugatedGradients.DecrementTolerance;
CGExitOptions.DecrementAbsTolerance=CtrlVar.Inverse.UaConjugatedGradients.DecrementAbsTolerance;
CGExitOptions.DecrementRelativeTolerance=CtrlVar.Inverse.UaConjugatedGradients.DecrementRelativeTolerance;
CGExitOptions.dJTolerance=CtrlVar.Inverse.UaConjugatedGradients.dJTolerance;

%% make initial iterate feasible
p=p(:);
p=kk_proj(p,pub,plb);

%%  Get the metric matrix
if CtrlVar.Inverse.RieszMapGradient
    G=MUA.G;
else
    G=1;
end

%% Evaluation cache, shared with the nested function PhiCached
CacheGamma=[] ; CacheJ=[] ; CacheGrad={} ; CacheFOuts={} ;
nFuncEval=0 ; nGradEval=0 ;

%% First evaluation of the cost function at the current point
[J0,dJdp,~,fOuts]=func(p);  nFuncEval=nFuncEval+1 ; nGradEval=nGradEval+1 ;
dJdp=dJdp(:);
mdJdp=-dJdp ; % this is the steepest descent direction in the G metric
d=mdJdp ;     % first search direction is simply the steepest-descent direction.

Gd=G*d ;
sGs=dJdp'*(G*dJdp);
Decrement = 0.5*sGs;
slope0=dJdp'*Gd;

GradNorm=sqrt(sGs); 
%% Make sure the RunInfo field is OK and properly set

if isempty(RunInfo) ||  ~isfield(RunInfo,'Inverse') || numel(RunInfo.Inverse.Iterations)<=1
    RunInfo.Inverse.Iterations(1)=0;
    RunInfo.Inverse.J(1)=J0;

    if CtrlVar.Inverse.StoreSolutionAtEachIteration
        RunInfo.Inverse.p{1}=p;
    end

    if isfield(fOuts,'R')
        RunInfo.Inverse.R(1)=fOuts.R;
    end
    if isfield(fOuts,'I')
        RunInfo.Inverse.I(1)=fOuts.I;
    end
    RunInfo.Inverse.StepSize(1)=0;
    RunInfo.Inverse.Decrement(1)=Decrement;
    RunInfo.Inverse.GradNorm(1)=GradNorm;
end

RunInfo.Inverse.ConjGradUpdate=0;

%% Is the starting point already stationary?
%
% If the gradient vanishes identically, as happens with perfect synthetic data and the true parameter field, then d=0
% and slope0=0, and the linear model below gives gamma1=-0.05*abs(J0)/slope0 = 0/0 = NaN. The first trial point would
% then be p+NaN*0 = NaN and the forward solve would fail. Catch it here instead.

if ~any(d) || ~(slope0<0)
    fprintf('\n +++++++++++ At the starting point the gradient is zero (decrement=%g, slope0=%g). \n',Decrement,slope0)
    fprintf(' +++++++++++ p is already a stationary point, nothing to do. \n\n')
    return
end

%% Determine an initial step size for the first line search
%
% A linear approximation aiming for a 5% reduction, refined by a quadratic through J0, slope0 and J at that trial step.
% The trial point is evaluated through PhiCached so that it is available to the first line search.
%
% Note the abs(J0): J is not guaranteed positive (JGH itself warns when J<0), and without it a negative J0 would give a
% negative trial step, i.e. a step in the ascent direction.

gamma1=-0.05*abs(J0)/slope0 ;

% J0 can be zero, or so small that the linear model provides no useful length scale, in which case gamma1 comes out as 0
% or NaN. Fall back on a step that moves p by roughly 0.1% of its own norm.
if ~isfinite(gamma1) || gamma1<=0
    gamma1=1e-3*max(norm(p),1)/max(norm(d),realmin) ;
end

[J1,~]=PhiCached(gamma1) ;

% the forward model may fail to converge for this step size, in which case reduce gamma until it does
while isnan(J1) && gamma1>eps
    fprintf("Objective function returned NaN; trying a new point...\n")
    gamma1=gamma1/10;
    [J1,~]=PhiCached(gamma1) ;
end

gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadratic approx

if ~isfinite(gamma) || gamma<=0  % quadratic estimate not useful, fall back on the linear one
    gamma=gamma1;
end

gammaStart=gamma;

%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g \t  decrement=%g \t \t gamma=%-g \n \n',J0,fOuts.I,fOuts.R,Decrement,gamma)

It0=RunInfo.Inverse.Iterations(end);

fprintf('\n   It\t #cgUpd F-count G-count   \t   J     \t   I    \t   R           decrement  \t gamma  \t SDratio    CGcorr \n')

fprintf('%5i\t%5i\t%5i\t%5i\t%15.10g\t%15.10g\t%15.10g\t  %15.10g \t %10g\t%8.4f %8.4f \n',...
    It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,nGradEval,J0,fOuts.I,fOuts.R,Decrement,gamma,nan,nan)

%%

ExitInfo=[];                              % state carried by CGExitCriteria

for Iteration=1:CtrlVar.Inverse.Iterations

    % J at gammaStart, together with its slope. The cache already holds this point whenever gammaStart was one of the
    % trial steps probed above, or by the previous line search.
    kc=find(CacheGamma==gammaStart,1,'last') ;
    if isempty(kc)
        [J1,slope1]=PhiCached(gammaStart) ;
    else
        J1=CacheJ(kc) ; slope1=CacheGrad{kc}'*Gd ;
    end

    % Supplying the slope at gammaStart lets LineSearchWolfe accept that step without any further evaluation whenever it
    % already satisfies both Wolfe conditions.
    LineSearchOptions.slopeStart=slope1 ;

    [gamma,~,LineSearchInfo]=LineSearchWolfe(slope0,gammaStart,J0,J1,@PhiCached,LineSearchOptions);

    gammaLastMinimum=gamma;

    pNew=p+gamma*d ;
    pProj=kk_proj(pNew,pub,plb);
    Projected=~isequal(pProj,pNew) ;
    p=pProj ;
    mdJdpLast=mdJdp;

    % The accepted point has already been evaluated inside the line search, so read it back rather than calling func
    % again. Only if kk_proj moved p, or the accepted step is somehow not in the cache, is a fresh evaluation needed.
    kc=find(CacheGamma==gamma,1,'last') ;
    if gamma==0
        % The line search found no point better than phi(0), so p is unchanged and J0, dJdp and fOuts still refer to it.
        % Nothing needs re-evaluating.
    elseif ~Projected && ~isempty(kc)
        J0=CacheJ(kc) ; dJdp=CacheGrad{kc} ; fOuts=CacheFOuts{kc} ;
    else
        [J0,dJdp,~,fOuts]=func(p) ;  nFuncEval=nFuncEval+1 ; nGradEval=nGradEval+1 ;
        dJdp=dJdp(:) ;
    end

    mdJdp=-dJdp ; % this is the steepest descent direction in the G metric

    sGs=dJdp'*(G*dJdp);
    Decrement = 0.5*sGs;
    GradNorm=sqrt(sGs);
    Misfit=fOuts.I;

  
    fprintf('%5i\t%5i\t%5i\t%5i\t%15.10g\t%15.10g\t%15.10g\t  %15.10g \t %10g\t%8.4f %8.4f \n',...
        Iteration+It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,nGradEval,J0,fOuts.I,fOuts.R,...
        Decrement,gamma,cgInfo.SufficientDescentRatio,cgInfo.CGCorrectionRatio)

    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.I];
    RunInfo.Inverse.Decrement=[RunInfo.Inverse.Decrement;Decrement];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];

    [Exit,ExitInfo]=CGExitCriteria(CGExitOptions,ExitInfo,Iteration,J0,sGs,cgInfo,LineSearchInfo,Misfit);

    if Exit
        fprintf('\n Inversion stopped after %i iterations, exit flag %i : %s \n',...
            Iteration,ExitInfo.Flag,ExitInfo.Message)
        break
    end

    % update search direction. ExitInfo.ForceRestart is set by CGExitCriteria when a line search has failed on a CG
    % direction, or on the first sign of stagnation, and asks for the CG history to be discarded.
    [d,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d,G,CtrlVar,cgInfo,ExitInfo.ForceRestart);

    RunInfo.Inverse.ConjGradUpdate=cgInfo.NumberOfConjGradUpdatesWithoutReset;

    slope0Last=slope0 ;
    Gd=G*d ;
    slope0=dJdp'*Gd;

    % Initial step for the next line search, Nocedal & Wright eq 3.60: scale the last accepted step by the ratio of the
    % directional derivatives at the two starting points. This is the standard choice for a Wolfe line search and is
    % better scaled than simply extending the previous step.
    if slope0<0 && isfinite(slope0Last) && slope0Last<0 && gammaLastMinimum>0
        gammaStart=gammaLastMinimum*slope0Last/slope0 ;
    else
        gammaStart=gammaLastMinimum ;
    end

    % gammaLastMinimum is ZERO whenever the line search found no point better than phi(0), i.e. LineSearchWolfe returned
    % iExit=-2. Falling back on it would feed gammaStart=0 into the next line search, which then errors on
    % ~(gammaStart>0). Re-derive a step from the linear model instead, which is positive because slope0<0.
    if ~isfinite(gammaStart) || gammaStart<=0
        gammaStart=-0.05*abs(J0)/slope0 ;
    end

    if ~isfinite(gammaStart) || gammaStart<=0
        fprintf('\n Inversion stopped after %i iterations: no positive step is available. \n',Iteration)
        break
    end

    % the cache refers to the old p and d, so it is no longer valid
    CacheGamma=[] ; CacheJ=[] ; CacheGrad={} ; CacheFOuts={} ;

end

% ExitInfo is still empty if CtrlVar.Inverse.Iterations<1, so the loop never ran
if isstruct(ExitInfo)
    fprintf('\n +++++++++++ Totals: %i cost function evaluations, %i gradient evaluations, %i CG restarts. \n\n',...
        nFuncEval,nGradEval,ExitInfo.nRestart)
else
    fprintf('\n +++++++++++ Totals: %i cost function evaluations, %i gradient evaluations. No iterations performed. \n\n',...
        nFuncEval,nGradEval)
end

%% nested functions

    function [J,slope]=PhiCached(gam)

        % Evaluates phi(gam)=J(p+gam*d) and its directional derivative, and caches everything the caller might need at
        % the accepted step. p, d and Gd are taken from the parent workspace and are fixed during a line search.

        [J,dJdpTrial,~,fOutsTrial]=func(p+gam*d) ;
        dJdpTrial=dJdpTrial(:) ;
        slope=dJdpTrial'*Gd ;

        nFuncEval=nFuncEval+1 ;
        nGradEval=nGradEval+1 ;

        CacheGamma(end+1,1)=gam ;
        CacheJ(end+1,1)=J ;
        CacheGrad{end+1,1}=dJdpTrial ;
        CacheFOuts{end+1,1}=fOutsTrial ;

    end

end


%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [d1,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo,ForceRestart)

if nargin<7 || isempty(ForceRestart) ; ForceRestart=false ; end

if ForceRestart
    % Discard the CG history and drop back onto steepest descent. This is the
    % same state the CG routine returns after one of its own resets.
    d1=mdJdp ;
    cgInfo.NumberOfConjGradUpdatesWithoutReset=0 ;
    cgInfo.teta=0 ;
    cgInfo.ConjGradAngle=0 ;
    cgInfo.ddAngle=nan ;
    cgInfo.SufficientDescentRatio=1 ;
    cgInfo.PowellRatio=nan ;
    cgInfo.CGCorrectionRatio=0 ;
    return
end

switch lower(CtrlVar.Inverse.UaConjugatedGradients.UpdateMethod)

    case {"-conjgrad-","conjgrad"}

        [d1,cgInfo]=NewConjugatedGradMetric(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo) ;

    case {"-steepestdecent-","steepestdecent"}

        d1=mdJdp ;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0 ;
        cgInfo.teta=0 ;
        cgInfo.ConjGradAngle=0 ;
        cgInfo.ddAngle=nan ;
        cgInfo.SufficientDescentRatio=1 ;
        cgInfo.PowellRatio=nan ;
        cgInfo.CGCorrectionRatio=0 ;

    otherwise

        error("CaseNotFound")

end

end
