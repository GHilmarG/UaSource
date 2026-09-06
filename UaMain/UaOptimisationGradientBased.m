function  [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%%
% This is basically a conjugated-gradient minimizer.
%
% It does a reasonably good job. Importantly it does allow for an arbitrary metric, which here is defined by the metric
% matrix G.  This is something that most optimization packages appear not to allow for.
%
% At the moment it is evaluating the cost function too often I think. I think the line-search, which is basically my modified
% backtracking algorithm with extrapolation allowed, needs improvements. Or maybe it should be replaced with a new line-search
% code that takes into account the Wolfe condition. The backtracking algorithm is really optimized for backtracking and only
% uses the Amarillo condition. 
%
% func is the function to me minimized
%  p is the parameter set, i.e. func(p)
%
%  Func is func evaluated as a function of step-size gamma in the direction of
%  the gradient: Func=@(gamma) func(p-gamma*dJdp);
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
JTolerance=0.01;
dJTolerance=0.0 ; 
dpTolerance=0.0;



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

GradNorm=norm(dJdp)/sqrt(numel(dJdp));

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
    RunInfo.Inverse.GradNorm=GradNorm;
    RunInfo.Inverse.ConjGradUpdate=0;
end


% determine initial search direction and initial step size for line-search.
if ~(isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0)
    gamma=CtrlVar.Inverse.InitialLineSearchStepSize;
    slope0=dJdp'*G*d;
else
    slope0=dJdp'*G*d;
    gamma1=-0.01*J0/slope0 ; % linear approx
    p1=p+gamma1*d;
    J1=func(p1);

    iCount=0 ; % sometimes the initial guess is so small that there is almost no change in the cost function
    % try to increase gamma by factor of 10 until at least 1% change has been generated.
    while (abs(J1-J0)/J1 < 0.01) && iCount<10

        gamma1=gamma1*10 ;
        p1=p+gamma1*d;
        J1=func(p1);
        iCount=iCount+1;
    end

    gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadratic approx
    if gamma<0 ; gamma=gamma1; end

end


%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \t \t gamma=%-g \n \n',J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,GradNorm,gamma)


Func=@(gamma) func(p+gamma*d);
J1=Func(gamma);
nFuncEval=1;


if isnan(J1)
    slope0=dJdp'*G*d; 
    gamma=-0.1*J0/slope0 ;  % modification on 23 Jan, 2019. Resetting gamma
    J1=Func(gamma);
    nFuncEval=nFuncEval+1;
end

while RunInfo.Forward.uvIterations==0

    % the gamma step caused so little change in the model parameters that the previous J0 uv solution was accepted.
    % So increase gamma
    fprintf(" Increasing the stepsize as the previous one caused insufficient changes in model parameters to require a new uv solution.\n")
    fprintf(" gamma increased from %g to %g \n",gamma,gamma*1000)
    gamma=gamma*1000 ;
    J1=Func(gamma);
    nFuncEval=nFuncEval+1;

end
gammaStart=gamma; 
%% Backtracking parameter modifications
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
%%

It0=RunInfo.Inverse.Iterations(end);

fprintf('\n   It   #fEval      J           I          R        decrement      gamma   #cgUpdate\n')
%fprintf('123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n')

fprintf('%5i  %5i %10g  %10g  %10g  %10g  %10g  %5i\n',It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,GradNorm,gamma,RunInfo.Inverse.ConjGradUpdate)
CtrlVar.GradientUpgradeMethod=CtrlVar.Inverse.GradientUpgradeMethod;

LineSearchInfo=[];
iBackTry=0; gammaLast=gamma;

for It=1:CtrlVar.Inverse.Iterations


    %CtrlVar.InfoLevelBackTrack=100 ;CtrlVar.doplots=1;
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    J1=Func(gammaStart);



    Gd = G*d ;                                   % once, outside the line search
    Phi = @(gamma) PhiEval(gamma,p,d,Gd,func) ;
    [gamma,JgammaNew,LineSearchInfo]=LineSearchWolfe(slope0,gammaStart,J0,J1,Phi,LineSearchInfo);
    nFuncEval=LineSearchInfo.nFuncEvaluations+1 ; % adding the one I did to get J1

    dJ=J0-JgammaNew; 
    % 
    % [gamma,JgammaNew,BackTrackingInfoVector]=BackTracking(slope0,gammaStart,J0,J1,Func,CtrlVar);
    % nFuncEval=nFuncEval+BackTrackingInfoVector.nFuncEval;


    dpNorm=gamma * sqrt(d'*G*d); 

    p=p+gamma*d;
    p=kk_proj(p,pub,plb);
    mdJdpLast=mdJdp;
 
    % Get the new steepest descent direction. 
    [J0,dJdp,~,fOuts]=func(p);
    mdJdp=-dJdp ; % this is the steepest descent direction in the G metric
    nFuncEval=nFuncEval+1; % here J0 and JgammaNew must be (almost) equal

 
    Decrement = sqrt(0.5*(dJdp'*G*dJdp)) ; % with dJdp the Riesz-mapped descent direction 

    if Decrement < CtrlVar.Inverse.DecrementTolerance

        fprintf("Breking out of inverse iteration loop as decrement has reached tolerance.\ n")
        fprintf("          Decrement=%g\n",Decrement)
        fprintf("Decrement tolerance=%g \n",CtrlVar.Inverse.DecrementTolerance)

        break 
    
    end

    
    if dJ<dJTolerance
        fprintf("dJ tolerance (%g) reached with dJ=%g. \n",dJTolerance,dJ)
        break
    end

    if J0<JTolerance  % I can use J0 here because it has already been updated, i.e. this is currently the best J
        fprintf("J tolerance (%g) reached with J=%g. \n",JTolerance,J)
        break
    end

    if dpNorm< dpTolerance
        fprintf("dp tolerance (%g) reached with dpNorm=%g. \n",dpTolerance,dpNorm)
        break
    end


    % update search direction.
    [d,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d,G,CtrlVar,cgInfo);
    
    RunInfo.Inverse.ConjGradUpdate=cgInfo.NumberOfConjGradUpdatesWithoutReset;

    Func=@(gamma) func(p+gamma*d);
    slope0=dJdp'*G*d;  

 

    gammaStart=2*gammaLast;


    fprintf('%5i  %5i %10g  %10g  %10g  %10g  \t %10g \t %5i\n',It+It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma,cgInfo.NumberOfConjGradUpdatesWithoutReset)

    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;Decrement];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];

end

end


%%%%%%%%%%%%%%%%



function [d1,cgInfo]=NextSearchDirection(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo)


switch lower(CtrlVar.GradientUpgradeMethod)

    case 'conjgrad'

        [d1,cgInfo]=NewConjugatedGradMetric(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo) ;

    otherwise

        d1=d0;
end

end



function [J,slope]=PhiEval(gamma,p0,d,Gd,func)
[J,dJdp]=func(p0+gamma*d) ;
slope=dJdp'*Gd ;
end


