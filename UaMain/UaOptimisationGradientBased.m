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
 
    CtrlVar.Inverse.UaConjugatedGradients.Armijo=1e-4;
    CtrlVar.Inverse.UaConjugatedGradients.WolfeCurvature=0.2;
    CtrlVar.Inverse.UaConjugatedGradients.MaxFuncEvalutionsInLineSearch=15;
    CtrlVar.Inverse.UaConjugatedGradients.InfoLevel=0;

    CtrlVar.Inverse.UaConjugatedGradients.UpdateMethod="-ConjGrad-" ; %{'SteepestDecent','ConjGrad'}
    CtrlVar.ConjugatedGradientsUpdate="FR" ; % {"FR","PR","HS",DY"}
    CtrlVar.Inverse.DecrementTolerance=1e-10;
end



%% Line-search parameters, here used for the LineSearchWolfe
LineSearchOptions.c1= CtrlVar.Inverse.UaConjugatedGradients.Armijo;  % Armijo parameter
LineSearchOptions.c2= CtrlVar.Inverse.UaConjugatedGradients.WolfeCurvature ;  % Wolfe, curvature parameter.
LineSearchOptions.MaxFuncEvaluations= CtrlVar.Inverse.UaConjugatedGradients.MaxFuncEvalutionsInLineSearch; 


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

if CtrlVar.Inverse.RieszMapGradient
    G=MUA.MetricMatrix; 
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

slope0=dJdp'*G*d;
gamma1=-0.05*J0/slope0 ; % linear approx, step size determined by hoping for 5% reduction based on slope at origin
p1=p+gamma1*d;
J1=func(p1);

% OK, so now I have J0 at gamma=0, J1 at gamma1, and the slope at gamma=0
% with these three numbers I can build a quadratic approximation and estimate the minimum gamma
gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadratic approx
if gamma<0 ; gamma=gamma1; end




%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  decrement=%g \t \t gamma=%-g \n \n',J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma)


J1=func(p+gamma*d);
nFuncEval=1;

while isnan(J1)
    gamma=gamma/10;
    J1=func(p+gamma*d);
    nFuncEval=nFuncEval+1;
end



%%
It0=RunInfo.Inverse.Iterations(end);
RunInfo.Inverse.ConjGradUpdatenFuncEval=0; 
nFuncEval=0;

%%
fprintf('\n   It #cgUpdate  F-count     J           I          R       decrement      gamma  \n')
%fprintf('123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n')


fprintf('%5i\t%5i\t%5i %10g  %10g  %10g  %10g  \t %10g \n',It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma)

%%

ExitInfo=[];                              % state carried by CGExitCriteria


gammaStart=gamma;

for Iteration=1:CtrlVar.Inverse.Iterations



    J1=func(p+gammaStart*d); nFuncEval=nFuncEval+1; 

    % This is the old backtracking approach. Now no longer used. Instead I'm using line search with Wolfe conditions
    %
    % CtrlVar.InfoLevelBackTrack=100 ; CtrlVar.doplots=1;   CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    % Func=@(gamma) func(p+gamma*d);
    % [gammaTest,JgammaTest,BackTrackingInfoVector]=BackTracking(slope0,gammaStart,J0,J1,Func,CtrlVar);
    %
    %
    %


    Gd = G*d ;                                   % once, outside the line search
    Phi = @(gamma) PhiEval(gamma,p,d,Gd,func) ;
    [gamma,JgammaNew,LineSearchInfo]=LineSearchWolfe(slope0,gammaStart,J0,J1,Phi,LineSearchOptions);
    nFuncEval=nFuncEval+LineSearchInfo.nFuncEvaluations; % adding the one I did to get J1

  
    gammaLastMinimum=gamma;
  
    p=p+gamma*d;
    p=kk_proj(p,pub,plb);
    mdJdpLast=mdJdp;

    % Get the new steepest descent direction.
    [J0,dJdp,~,fOuts]=func(p);   nFuncEval=nFuncEval+1; % 
    mdJdp=-dJdp ; % this is the steepest descent direction in the G metric
    nFuncEval=nFuncEval+1; % 


    sGs=dJdp'*G*dJdp;
    Decrement = 0.5*sGs;
    Misfit=fOuts.MisfitOuts.I;

    % Record this iteration BEFORE testing for exit, otherwise the last
    % iteration is missing from RunInfo whenever the loop breaks.
   
    fprintf('%5i\t%5i\t%5i %10g  %10g  %10g  %10g  \t %10g \n',Iteration+It0,cgInfo.NumberOfConjGradUpdatesWithoutReset,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,Decrement,gamma)

    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;Decrement];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];
    
    CtrlVar.Inverse.InfoLevel=0; 
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

switch lower(CtrlVar.Inverse.UaConjugatedGradients.UpdateMethod)

    case {"-conjgrad-","conjgrad"}

        [d1,cgInfo]=NewConjugatedGradMetric(mdJdp,mdJdpLast,d0,G,CtrlVar,cgInfo) ;

    case {"-steepestdecent-","steepestdecent"}

    
        d1=mdJdp ;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0 ;
        cgInfo.teta=0 ;
        cgInfo.ConjGradAngle=0 ;
        cgInfo.ddAngle=nan ;

    otherwise

        error("CaseNotFound")

end

end



function [J,slope]=PhiEval(gamma,p0,d,Gd,func)
[J,dJdp]=func(p0+gamma*d) ;
slope=dJdp'*Gd ;
end


