function  [p,UserVar,RunInfo]=UaOptimisationGradientBasedBounded(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%%
% Bounded version of UaOptimisationGradientBased.
%
% Identical to that routine except that the box constraints plb <= p <= pub are respected by the line search itself,
% rather than by projecting the accepted step afterwards.
%
%
%% Why this is needed
%
% Projecting the arc, i.e. minimising phi(gamma)=J(P(p+gamma*d)), gives a piecewise-smooth function with a kink wherever
% a component of p reaches a bound. A Wolfe line search assumes phi is C^1, so at a kink the curvature condition is
% meaningless, and conjugacy is destroyed because the direction effectively changes subspace in the middle of a step.
%
% Instead, the search direction is restricted to the current face of the box, and the step is truncated at the first
% bound that would be crossed:
%
%   1) components of d pointing out of the box at a variable already sitting on a bound are set to zero. Without this,
%      gammaMax would be zero and no step would be possible at all.
%
%   2) gammaMax is the largest gamma for which p+gamma*d stays feasible, and is passed to LineSearchWolfe.
%
% phi is then smooth on [0,gammaMax] and the Wolfe conditions mean what they normally mean. LineSearchWolfe returning
% iExit=4 now signals that a new bound has become active, which is information rather than failure.
%
%   3) whenever the active set changes, the CG history is discarded. This restart is mandatory, not a safeguard: the
%      accumulated conjugacy refers to a subspace that no longer exists.
%
%   4) the decrement is formed from the PROJECTED gradient. The full gradient does not vanish at a constrained solution,
%      since the active components are held by the bounds, so a decrement built from it would never reach its tolerance
%      and every bounded inversion would run to the iteration limit.
%
%
%% The free-variable descent direction
%
% The steepest-descent direction on the current face is NOT obtained by masking the Riesz-mapped gradient. With
% s = -M.*(G\g), where M is the free-variable mask, the slope is
%
%   g's = -(M.*r)'G(M.*r) - ((1-M).*r)'G(M.*r) ,   r = G\g
%
% and the second term has no definite sign, so descent fails once many variables are active. That was verified: on a
% bound-constrained quadratic with 48 of 60 variables active at the solution, masking produced an ascent direction after
% 44 iterations and the run had to abort.
%
% Instead the Riesz map is restricted to the free variables, i.e. G_FF s_F = -g_F is solved directly. Then
%
%   g's = -g_F' (G_FF \ g_F) < 0
%
% so descent is guaranteed. Because s vanishes on the active set, every inner product x'Gy formed inside
% NewConjugatedGradMetric automatically reduces to x_F' G_FF y_F, so the CG update needs no modification at all. The
% factorisation of G_FF is cached and rebuilt only when the active set changes.
%
%
%% The working set
%
% A variable enters the active set as soon as it sits on a bound and the gradient pushes it further out. It is released
% only once the gradient pushes it back INTO the box by more than a threshold gRelease. Between those two conditions,
% i.e. for -gRelease < g < 0 at a lower bound, the variable keeps whatever state it had.
%
% That hysteresis is what stops the active set oscillating. Testing the sign of g alone makes a node whose gradient
% hovers near zero flip state on alternate iterations, and because every set change forces a CG restart, the conjugate
% gradients never accumulate. This was measured on a Hoffsjokull B inversion: over the last 200 of 1000 iterations, 95%
% of steps changed nActive by exactly one and 92% of those alternated in sign, with gamma and the decrement oscillating
% with period two. Only 1.6% of iterations in that block managed a CG update, and 481 entries against 402 exits produced
% a net change of just 82.
%
% Position alone cannot fix this: after truncation the variable sits exactly on its bound, so any position tolerance
% classifies it identically on both iterations. The oscillation is in g, so the remedy has to be in g.
%
% gRelease is scaled by the root-mean-square of the free-variable gradient, so it shrinks as the iteration converges.
% A fixed threshold would eventually freeze variables that genuinely want to leave the bound.
%
% The test uses the l2 gradient g=G*dJdp, since the KKT conditions are stated in terms of g and not of its Riesz
% representer.
%
%
%% Cost-function evaluations
%
% As in UaOptimisationGradientBased: every evaluation goes through the nested function PhiCached, which stores gamma, J,
% dJdp and the JGH outputs for each trial step, so the accepted point need not be re-evaluated. Because the step is now
% feasible by construction, kk_proj almost never alters p, so the cache is used essentially every iteration.
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

    CtrlVar.Inverse.UaConjugatedGradients.ActiveSetPositionTolerance=1e-8;
    CtrlVar.Inverse.UaConjugatedGradients.ActiveSetReleaseTolerance=1e-2;
end

% Active-set tolerances. See the note on the working set in the header.
%
%   ActiveSetPositionTolerance : how close to a bound counts as being on it, relative to the box width
%   ActiveSetReleaseTolerance  : how strongly the gradient must push back into the box before a variable is released,
%                                relative to the rms free-variable gradient. Set to 0 to recover the old behaviour of
%                                testing the sign of g alone, which oscillates.

epsPos=GetOpt(CtrlVar.Inverse.UaConjugatedGradients,'ActiveSetPositionTolerance',1e-8) ;
epsRelease=GetOpt(CtrlVar.Inverse.UaConjugatedGradients,'ActiveSetReleaseTolerance',1e-2) ;

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
plb=plb(:) ; pub=pub(:) ;
p=kk_proj(p,pub,plb);

% Tolerance used to decide whether a variable sits on a bound. Relative, because box widths can span several orders of
% magnitude: for B the width is about 9.3*(s-S), which varies from metres at the ice margin to many kilometres inland.
BoxTol=epsPos*max(1,pub-plb) ;

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

GFFdec=[] ; GFFmask=[] ;                       % cached factorisation of G(Free,Free)

g=Gmul(G,dJdp) ;                               % l2 gradient
ActiveMask=BoxWorkingSet(p,g,plb,pub,BoxTol) ; % first call: no previous state, so no hysteresis yet
nEnter=0 ; nExit=0 ;
[mdJdp,sGs,GFFdec,GFFmask]=FreeSpaceDescent(g,ActiveMask,G,GFFdec,GFFmask) ;

d=mdJdp ;     % first search direction is the steepest-descent direction on the current face
[d,gammaMax]=BoundTruncate(p,d,plb,pub,BoxTol) ;

Gd=Gmul(G,d) ;
Decrement = 0.5*sGs;
slope0=g'*d;

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
% With perfect synthetic data, and B set to the true field, both I and R are zero and the gradient vanishes identically.
% The projected gradient is then zero, the search direction d is zero, and gammaMax comes out as Inf because there is no
% component of d of either sign. The linear model for the first step then gives
%
%       gamma1 = -0.05*abs(J0)/slope0 = 0/0 = NaN
%
% and since MATLAB's min ignores NaN, min(NaN,0.99*Inf) returns Inf rather than NaN. The first trial point is then
% p+Inf*0 = NaN and the forward solve fails. Catch the stationary point here instead.

if ~any(d) || ~(slope0<0)
    fprintf('\n +++++++++++ At the starting point the projected gradient is zero (decrement=%g, slope0=%g). \n',Decrement,slope0)
    fprintf(' +++++++++++ p is already a stationary point of the bound-constrained problem, nothing to do. \n\n')
    return
end

%% Determine an initial step size for the first line search

gamma1=-0.05*abs(J0)/slope0 ;

% J0 can be zero, or so small that the linear model provides no useful length scale, in which case gamma1 comes out as 0
% or NaN. Fall back on a step that moves p by roughly 0.1% of its own norm. This guard MUST come before the min() below:
% min(NaN,x) returns x, so a NaN would otherwise vanish silently rather than being caught.
if ~isfinite(gamma1) || gamma1<=0
    gamma1=1e-3*max(norm(p),1)/max(norm(d),realmin) ;
end

gamma1=min(gamma1,0.99*gammaMax) ;          % keep the first probe strictly inside the box
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

gammaStart=min(gamma,0.99*gammaMax);

%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  decrement=%g \t \t gamma=%-g \t active=%i of %i \n \n',...
    J0,fOuts.I,fOuts.R,Decrement,gammaStart,sum(ActiveMask),numel(p))

It0=RunInfo.Inverse.Iterations(end);

fprintf('\n   It\t #cgUpd F-count G-count   \t   J     \t   I    \t   R  \t  decrement        gamma  \t SDratio    CGcorr  nActive  nFlip \n')

fprintf('%5i\t%5i\t%5i\t%5i\t%15.10g\t%15.10g\t%15.10g\t%10g \t %10g\t%8.4f %8.4f %8i %6i \n',...
    It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,nGradEval,J0,fOuts.I,fOuts.R,...
    Decrement,gammaStart,nan,nan,sum(ActiveMask),0)

%%

ExitInfo=[];                              % state carried by CGExitCriteria

for Iteration=1:CtrlVar.Inverse.Iterations

    % J at gammaStart, together with its slope.
    kc=find(CacheGamma==gammaStart,1,'last') ;
    if isempty(kc)
        [J1,slope1]=PhiCached(gammaStart) ;
    else
        J1=CacheJ(kc) ; slope1=CacheGrad{kc}'*Gd ;
    end

    LineSearchOptions.slopeStart=slope1 ;
    LineSearchOptions.gammaMax=gammaMax ;     % <- the step is never allowed out of the box

    [gamma,~,LineSearchInfo]=LineSearchWolfe(slope0,gammaStart,J0,J1,@PhiCached,LineSearchOptions);

    gammaLastMinimum=gamma;

    pNew=p+gamma*d ;
    pProj=kk_proj(pNew,pub,plb);
    % gamma<=gammaMax makes pNew feasible by construction, so kk_proj should only ever clip at roundoff level. A strict
    % isequal test would then flag a spurious projection and throw away the cached evaluation.
    Projected = max(abs(pProj-pNew)) > 1e-12*max(1,max(abs(pNew))) ;
    p=pProj ;
    mdJdpLast=mdJdp;

    kc=find(CacheGamma==gamma,1,'last') ;
    if gamma==0
        % The line search found no point better than phi(0), so p is unchanged and J0, dJdp and fOuts still refer to it.
        % Nothing needs re-evaluating. CGExitCriteria sees LineSearchWolfe iExit=-2 and will either request a restart on
        % steepest descent, or stop if the direction already was steepest descent.
    elseif ~Projected && ~isempty(kc)
        J0=CacheJ(kc) ; dJdp=CacheGrad{kc} ; fOuts=CacheFOuts{kc} ;
    else
        [J0,dJdp,~,fOuts]=func(p) ;  nFuncEval=nFuncEval+1 ; nGradEval=nGradEval+1 ;
        dJdp=dJdp(:) ;
    end

    % Working set at the new point, and the steepest-descent direction on that face.
    g=Gmul(G,dJdp) ;
    ActiveMaskLast=ActiveMask ;
    gRelease=ReleaseThreshold(g,ActiveMaskLast,epsRelease) ;
    ActiveMask=BoxWorkingSet(p,g,plb,pub,BoxTol,ActiveMaskLast,gRelease) ;
    SetChanged=~isequal(ActiveMask,ActiveMaskLast) ;

    nEnter=nEnter+sum(ActiveMask & ~ActiveMaskLast) ;
    nExit =nExit +sum(~ActiveMask &  ActiveMaskLast) ;
    nFlip =sum(ActiveMask~=ActiveMaskLast) ;

    [mdJdp,sGs,GFFdec,GFFmask]=FreeSpaceDescent(g,ActiveMask,G,GFFdec,GFFmask) ;

    Decrement = 0.5*sGs;
    GradNorm=sqrt(sGs);
    Misfit=fOuts.I;

    fprintf('%5i\t%5i\t%5i\t%5i\t%15.10g\t%15.10g\t%15.10g\t%10g \t %10g\t%8.4f %8.4f %8i %6i \n',...
        Iteration+It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,nGradEval,J0,fOuts.I,fOuts.R,...
        Decrement,gamma,cgInfo.SufficientDescentRatio,cgInfo.CGCorrectionRatio,sum(ActiveMask),nFlip)

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

    % Restart the CG whenever the active set has changed, in addition to the restarts requested by CGExitCriteria.
    % iExit=4 means the line search ran into gammaMax, i.e. a new bound was reached.
    ForceRestart = ExitInfo.ForceRestart || SetChanged || LineSearchInfo.iExit==4 ;

    [d,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d,G,CtrlVar,cgInfo,ForceRestart);

    RunInfo.Inverse.ConjGradUpdate=cgInfo.NumberOfConjGradUpdatesWithoutReset;

    [d,gammaMax]=BoundTruncate(p,d,plb,pub,BoxTol) ;

    slope0Last=slope0 ;
    Gd=Gmul(G,d) ;
    slope0=g'*d;

    % Descent guard. With the free-variable Riesz map this should not trigger, since g'*mdJdp<0 is guaranteed and
    % NewConjugatedGradMetric enforces sufficient descent in the same metric. Kept as a safety net.
    if ~(slope0<0) && any(d)
        [d,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d,G,CtrlVar,cgInfo,true);
        [d,gammaMax]=BoundTruncate(p,d,plb,pub,BoxTol) ;
        Gd=Gmul(G,d) ;
        slope0=g'*d;
        if CtrlVar.Inverse.UaConjugatedGradients.InfoLevel>=1
            fprintf('   CG direction was not a descent direction, restarted on steepest descent. \n')
        end
    end

    if ~any(d)
        fprintf('\n Inversion stopped after %i iterations: the projected gradient is zero, all variables are at a bound. \n',Iteration)
        break
    end

    if ~(slope0<0)
        fprintf('\n Inversion stopped after %i iterations: no descent direction could be found on the current face (slope0=%g). \n',Iteration,slope0)
        break
    end

    % Initial step for the next line search, Nocedal & Wright eq 3.60.
    if slope0<0 && isfinite(slope0Last) && slope0Last<0 && gammaLastMinimum>0
        gammaStart=gammaLastMinimum*slope0Last/slope0 ;
    else
        gammaStart=gammaLastMinimum ;
    end

    % gammaLastMinimum is ZERO whenever the line search found no point better than phi(0), i.e. LineSearchWolfe returned
    % iExit=-2. Falling back on it would feed gammaStart=0 into the next line search, which then errors on
    % ~(gammaStart>0). Re-derive a step from the linear model instead: this is positive because slope0<0 has just been
    % checked above. The same applies if the ratio in eq 3.60 came out non-finite.
    if ~isfinite(gammaStart) || gammaStart<=0
        gammaStart=-0.05*abs(J0)/slope0 ;
    end

    gammaStart=min(gammaStart,0.99*gammaMax) ;   % the first probe must stay inside the box

    if ~isfinite(gammaStart) || gammaStart<=0
        fprintf('\n Inversion stopped after %i iterations: no positive feasible step is available (gammaMax=%g). \n',...
            Iteration,gammaMax)
        break
    end

    % the cache refers to the old p and d, so it is no longer valid
    CacheGamma=[] ; CacheJ=[] ; CacheGrad={} ; CacheFOuts={} ;

end

if isstruct(ExitInfo)
    fprintf('\n +++++++++++ Totals: %i cost function evaluations, %i gradient evaluations, %i CG restarts, %i variables at a bound. \n',...
        nFuncEval,nGradEval,ExitInfo.nRestart,sum(ActiveMask))
    fprintf(' +++++++++++ Active set: %i entries, %i exits, net %+i. Churn (entries+exits-|net|) = %i. \n\n',...
        nEnter,nExit,nEnter-nExit,nEnter+nExit-abs(nEnter-nExit))
else
    fprintf('\n +++++++++++ Totals: %i cost function evaluations, %i gradient evaluations. No iterations performed. \n\n',...
        nFuncEval,nGradEval)
end

%% nested functions

    function [J,slope]=PhiCached(gam)

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


function [d,gammaMax,nZeroed]=BoundTruncate(p,d,plb,pub,BoxTol)

% Restrict the search direction to the current face of the box, and return the largest step for which p+gamma*d remains
% feasible.
%
% Components sitting on a bound with d pointing out of the box are set to zero. Without this gammaMax would be zero and
% no step would be possible. Note that variables with a degenerate box, plb=pub, are frozen automatically: they are both
% AtLower and AtUpper, so d is zeroed there whatever its sign.

AtLower = p <= plb+BoxTol ;
AtUpper = p >= pub-BoxTol ;

Out = (AtLower & d<0) | (AtUpper & d>0) ;
nZeroed = sum(Out) ;
d(Out)=0 ;

gammaMax=inf ;

i = d>0 ;
if any(i) ; gammaMax=min(gammaMax,min((pub(i)-p(i))./d(i))) ; end

i = d<0 ;
if any(i) ; gammaMax=min(gammaMax,min((plb(i)-p(i))./d(i))) ; end

gammaMax=max(gammaMax,0) ;   % guard against a small negative value from roundoff

end

%%

function Active=BoxWorkingSet(p,g,plb,pub,BoxTol,ActiveLast,gRelease)

% The working set, with hysteresis. See the note in the header.
%
%   enter   : sitting on a bound and the gradient pushes further out
%   release : sitting on a bound and the gradient pushes back in by more than gRelease
%   retain  : in between, keep the previous state
%
% Testing the sign of g alone (gRelease=0) makes a variable whose gradient hovers near zero flip state on alternate
% iterations, which forces a CG restart every iteration and prevents the conjugate gradients from ever accumulating.
%
% Variables with a degenerate box, plb=pub, satisfy both position tests and are therefore always active, i.e. frozen.

if nargin<6 || isempty(ActiveLast) ; ActiveLast=false(size(p)) ; end
if nargin<7 || isempty(gRelease)   ; gRelease=0 ; end

AtLower = p <= plb+BoxTol ;
AtUpper = p >= pub-BoxTol ;
AtBound = AtLower | AtUpper ;

PushOut      = (AtLower & g >  0       ) | (AtUpper & g < 0        ) ;
PushInStrong = (AtLower & g < -gRelease) | (AtUpper & g >  gRelease) ;

Active = AtBound & ( PushOut | (ActiveLast & ~PushInStrong) ) ;

Active = Active | (AtLower & AtUpper) ;     % degenerate box, frozen

end

%%

function gRelease=ReleaseThreshold(g,ActiveMask,epsRelease)

% Threshold the gradient must exceed, pointing back into the box, before a variable is released from the active set.
%
% Scaled by the rms of the free-variable gradient so that it shrinks as the iteration converges: a fixed threshold would
% eventually freeze variables that genuinely want to leave their bound.

if epsRelease<=0 ; gRelease=0 ; return ; end

Free=~ActiveMask ;
if any(Free) ; gScale=sqrt(mean(g(Free).^2)) ; else ; gScale=sqrt(mean(g.^2)) ; end

gRelease=epsRelease*gScale ;

end

%%

function v=GetOpt(S,Field,Default)

% value of S.(Field) if it is present and non-empty, otherwise Default

if isstruct(S) && isfield(S,Field) && ~isempty(S.(Field))
    v=S.(Field) ;
else
    v=Default ;
end

end

%%

function y=Gmul(G,x)

% G may be the metric matrix, or the scalar 1 when CtrlVar.Inverse.RieszMapGradient is false. Both cases are handled by
% the same expression; this exists only so that the intent is explicit at the call sites.

y=G*x ;

end

%%

function [s,sGs,GFFdec,GFFmask]=FreeSpaceDescent(g,ActiveMask,G,GFFdec,GFFmask)

% Steepest-descent direction on the current face, in the G metric, together with the corresponding decrement.
%
% Solves  G_FF s_F = -g_F  and sets s=0 on the active set, so that
%
%       g's   = -g_F' (G_FF \ g_F) = -sGs < 0
%
% which guarantees descent, and
%
%       s'G s = g_F' (G_FF \ g_F) = sGs
%
% so sGs is both the squared dual norm of the projected gradient and twice the decrement. It vanishes exactly at a KKT
% point of the bound-constrained problem, which is what makes the convergence test in CGExitCriteria meaningful here.
%
% The factorisation of G_FF is cached and rebuilt only when the active set changes.

s=zeros(size(g)) ;
sGs=0 ;

Free=~ActiveMask ;
if ~any(Free) ; return ; end

if isscalar(G)

    s(Free)=-g(Free)/G ;

else

    if isempty(GFFmask) || ~isequal(GFFmask,ActiveMask)
        GFF=G(Free,Free) ;
        GFF=(GFF+GFF')/2 ;      % a principal submatrix of an SPD matrix is SPD; symmetrise against roundoff
        GFFdec=decomposition(GFF,'chol') ;
        GFFmask=ActiveMask ;
    end

    s(Free)=-(GFFdec\g(Free)) ;

end

sGs=max(-g(Free)'*s(Free),0) ;

end

%%

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
