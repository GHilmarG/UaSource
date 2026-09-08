function  [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%%
% This is basically a conjugated-gradient minimizer.
%
% It does a reasonably good job. Importantly it does allow for an arbitrary metric, which here is defined by the metric
% matrix G.  This is something that most optimization packages appear not to allow for.
%
% The line search is now done using a new function,LineSearchWolfe, and it returns a minimum satisfying both Wolfe
% conditions, i.e. both the Armijo rule and the curvature condition. 
% 
%
%
% func is the function to me minimized
%
%  p is the parameter set, i.e. func(p)
%
%
%%



narginchk(8,8)
nargoutchk(3,3)

if isempty(CtrlVar)

    CtrlVar.Inverse.InitialLineSearchStepSize=1;
    CtrlVar.Inverse.HessianEstimate='0';
    CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize=1e-5;
    CtrlVar.Inverse.MinimumRelativelLineSearchStepSize=1e-4;
    CtrlVar.Inverse.MaximumNumberOfLineSeachSteps=100;
    CtrlVar.Inverse.InfoLevel=1;
    CtrlVar.Inverse.Iterations=4;
    CtrlVar.Inverse.GradientUpgradeMethod="-ConjGrad-" ; %{'SteepestDecent','ConjGrad'}
    CtrlVar.NewtonAcceptRatio=0.5;
    CtrlVar.NLtol=1e-15;
    CtrlVar.doplots=1;
    CtrlVar.Inverse.StoreSolutionAtEachIteration=0;
end

CtrlVar.Inverse.DecrementTolerance=1e-10;



cgInfo.NumberOfConjGradUpdatesWithoutReset=0;
cgInfo.ConjGradAngle=nan;
cgInfo.ddAngle=nan;
%%
p=p(:);
p=kk_proj(p,pub,plb);


%% Note: for anything other then the l2 gradient, this is not quite OK
%
% To do: I need to make sure that I'm using the same metric matrix, G, throughout the code.
% This is now easy to do as I calculate the metric matrix in one function, but I still need to change the
% ApplyAdjointGradientPreMultiplier.m
%
% I need the MetricMatrix, for example for L2 this would be
if CtrlVar.Inverse.AdjointGradientPreMultiplier=="M" || CtrlVar.Inverse.AdjointGradientPreMultiplier=="L2"

    G=MUA.M/MUA.Area ;
    [isA,isB,isC] = isABC(CtrlVar);

    if isA && isC
        G=blkdiag(G,G);
    end
elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="H1"
    error("not implemented")
else
    G=1;
end


[J0,dJdp,~,fOuts]=func(p);
dJdp=dJdp(:);
mdJdp=-dJdp ; % this is the steepest descent direction in the G metric
d=mdJdp ;     % first search direction is simply the steepest-descent direction.

sGs=dJdp'*G*dJdp;
Decrement = 0.5*sGs;

RunInfo.Inverse.ConjGradUpdate=0;

if isempty(RunInfo) ||  numel(RunInfo.Inverse.Iterations)<=1
    RunInfo.Inverse.Iterations(1)=0;
    RunInfo.Inverse.J(1)=J0;

    if CtrlVar.Inverse.StoreSolutionAtEachIteration
        RunInfo.Inverse.p{1}=p;
    end

    if isfield(fOuts,'R')
        RunInfo.Inverse.R(1)=fOuts.RegOuts.R;
    end
    if isfield(fOuts,'I')
        RunInfo.Inverse.I(1)=fOuts.MisfitOuts.I;
    end
    RunInfo.Inverse.StepSize(1)=0;
    RunInfo.Inverse.GradNorm=Decrement;
    RunInfo.Inverse.ConjGradUpdate=0;
end


% determine initial search direction and initial step size for line-search.
if ~(isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0)
    gamma=CtrlVar.Inverse.InitialLineSearchStepSize;
    slope0=dJdp'*G*d;
else
    slope0=dJdp'*G*d;
    gamma1=-0.05*J0/slope0 ; % linear approx
    p1=p+gamma1*d;
    J1=func(p1);

    % OK, so now I have J0 at gamma=0, J1 at gamma1, and the slope at gamma=0
    % with these three numbers I can build a quadratic approximation and estimate the minimum gamma
    gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadratic approx
    if gamma<0 ; gamma=gamma1; end

end


%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \t \t gamma=%-g \n \n',J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma)


J1=func(p+gamma*d);
nFuncEval=1;

while isnan(J1)
    gamma=gamma/10;
    J1=func(p+gamma*d);
    nFuncEval=nFuncEval+1;
end


%% Backtracking parameter modifications
% These are here because I used to use the backtracking algorithm, but now I'm using the line-search algorithm, so these are
% now reduntant and can be deleted in the future.
%
CtrlVar.BacktrackingGammaMin=CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize;
CtrlVar.BackTrackMinXfrac=CtrlVar.Inverse.MinimumRelativelLineSearchStepSize;
CtrlVar.BackTrackMaxIterations=CtrlVar.Inverse.MaximumNumberOfLineSearchSteps;
CtrlVar.InfoLevelBackTrack=CtrlVar.Inverse.InfoLevelBackTrack;
CtrlVar.BackTrackGuardLower=0.25;
CtrlVar.BackTrackGuardUpper=0.95;
% Backtracking continues even if target has been reached if last reduction in
% ratio is smaller than:
CtrlVar.BackTrackContinueIfLastReductionRatioLessThan=0.5;
CtrlVar.NewtonAcceptRatio=0.5;
CtrlVar.BackTrackExtrapolationRatio=10;


%% Line-search parameters, here used for the LineSearchWolfe
LineSearchOptions.c1=1e-4;  % Armijo parameter
LineSearchOptions.c2=0.2 ;  % Wolfe, curvature parameter.



%%
It0=RunInfo.Inverse.Iterations(end);

fprintf('\n   It   #fEval      J           I          R        decrement      gamma  \t #cgUpdate\n')
%fprintf('123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n')

fprintf('%5i  %5i %10g  %10g  %10g  %10g  %10g  %5i\n',It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma,RunInfo.Inverse.ConjGradUpdate)
CtrlVar.GradientUpgradeMethod=CtrlVar.Inverse.GradientUpgradeMethod;


ExitInfo=[];                              % state carried by CGExitCriteria


gammaStart=gamma;

for Iteration=1:CtrlVar.Inverse.Iterations


    %CtrlVar.InfoLevelBackTrack=100 ;CtrlVar.doplots=1;
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;

    J1=func(p+gammaStart*d);



    Gd = G*d ;                                   % once, outside the line search
    Phi = @(gamma) PhiEval(gamma,p,d,Gd,func) ;
    [gamma,JgammaNew,LineSearchInfo]=LineSearchWolfe(slope0,gammaStart,J0,J1,Phi,LineSearchOptions);
    nFuncEval=LineSearchInfo.nFuncEvaluations+1 ; % adding the one I did to get J1

  
    gammaLastMinimum=gamma;
    %
    % [gamma,JgammaNew,BackTrackingInfoVector]=BackTracking(slope0,gammaStart,J0,J1,Func,CtrlVar);
    % nFuncEval=nFuncEval+BackTrackingInfoVector.nFuncEval;

    p=p+gamma*d;
    p=kk_proj(p,pub,plb);
    mdJdpLast=mdJdp;

    % Get the new steepest descent direction.
    [J0,dJdp,~,fOuts]=func(p);  % this is with gamma=0
    mdJdp=-dJdp ; % this is the steepest descent direction in the G metric
    nFuncEval=nFuncEval+1; % here J0 and JgammaNew must be (almost) equal


    sGs=dJdp'*G*dJdp;
    Decrement = 0.5*sGs;
    Misfit=fOuts.MisfitOuts.I;

    % Record this iteration BEFORE testing for exit, otherwise the last
    % iteration is missing from RunInfo whenever the loop breaks.
    fprintf('%5i  %5i %10g  %10g  %10g  %10g  \t %10g \t %5i\n',Iteration+It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma,cgInfo.NumberOfConjGradUpdatesWithoutReset)

    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;Decrement];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];

    [Exit,ExitInfo]=CGExitCriteria(CtrlVar,ExitInfo,Iteration,J0,sGs,cgInfo,LineSearchInfo,Misfit);

    if Exit
        fprintf('\n Inversion stopped after %i iterations, exit flag %i : %s \n',...
            Iteration,ExitInfo.Flag,ExitInfo.Message)
        fprintf(' Totals: %i cost function evaluations, %i gradient evaluations, %i CG restarts. \n\n',...
            ExitInfo.nFuncTotal,ExitInfo.nGradTotal,ExitInfo.nRestart)
        break
    end

  

    % update search direction. ExitInfo.ForceRestart is set by CGExitCriteria
    % when a line search has failed on a CG direction, or on the first sign of
    % stagnation, and asks for the CG history to be discarded.
    [d,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d,G,CtrlVar,cgInfo,ExitInfo.ForceRestart);

    RunInfo.Inverse.ConjGradUpdate=cgInfo.NumberOfConjGradUpdatesWithoutReset;

    slope0=dJdp'*G*d;

    % What is a sensible start value for gamma for the next line-search?
    % Should I extend it a bit from last minimum, or since I'm now using a line search maybe best to just use directly
    % gammaLastMinimum? 
    gammaStart=2*gammaLastMinimum;

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
    return
end

switch lower(CtrlVar.GradientUpgradeMethod)

    case 'conjgrad'

        [d1,cgInfo]=NewConjugatedGradMetric(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo) ;

    otherwise

        % Steepest descent. Note that this previously returned d1=d0, i.e. the
        % PREVIOUS search direction, which froze the direction at the initial
        % -dJdp for the whole run, so that every subsequent line search was a
        % search along a stale direction.
        d1=mdJdp ;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0 ;
        cgInfo.teta=0 ;
        cgInfo.ConjGradAngle=0 ;
        cgInfo.ddAngle=nan ;
end

end



function [J,slope]=PhiEval(gamma,p0,d,Gd,func)
[J,dJdp]=func(p0+gamma*d) ;
slope=dJdp'*Gd ;
end


