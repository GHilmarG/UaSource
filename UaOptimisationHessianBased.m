function  [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)
%
% func is the function to me minimized
%  p is the paramter set, i.e. func(p)
%
%  Func is func evaluated as a function of stepsize gamma in the direction of
%  the gradient: Func=@(gamma) func(p-gamma*dJdp);
%

narginchk(8,8)
nargoutchk(3,3)

%%
p=p(:);
p=kk_proj(p,pub,plb);

[J,dJdp,Hess,fOuts]=func(p); RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;


GradNorm=norm(dJdp);
RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
RunInfo.Inverse.J=[RunInfo.Inverse.J;J];
RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;NaN];




fprintf("   It \t #fEval\t    gamma   \t\t    J \t\t\t\t  J/J0 \t\t\t   |dp|/|p| \t\t |dJ/dp| \t sub-obtimality gap\n")
fprintf("%5i \t %5i \t %7.4g  \t\t %10.5g \t\t  %10.5g \t\t %10.5g \t\t %10.5g \t\t %g \n",0,RunInfo.Inverse.nFuncEval,NaN,J,NaN,NaN,GradNorm,NaN)

NewtonAcceptRatioMin=0.9 ;
CtrlVar.NewtonAcceptRatio=NewtonAcceptRatioMin ;
CtrlVar.BackTrackMinXfrac=1e-3 ; 
gamma=NaN;
for Iteration=1:CtrlVar.Inverse.Iterations
    
    J0=J;
    
    % [J0,dJdp,Hess]=func(p);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
    
    dp=-Hess\dJdp ;   % Newton system
    
    SubOptimality=-dJdp'*dp/2  ;  % sub-optimality at the beginning of the iteration step
    slope0=dJdp'*dp;
    
    
    % A sensible step length is something like
    if isnan(gamma)
        gamma=-0.01*J0/slope0  ;
    end
    
    theta=acos(-dJdp'*dp/(norm(dJdp)*norm(dp))); 
     fprintf("angle between search direction and steepest-decent direction is %f degrees.\n",theta*180/pi)
    
    [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
    
    % Test slope
    
    
    Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp
    
    TestSlope(Func,0,0.1*gamma,slope0) ;
    
    
    
    
    
    %slope0=[];
    
    BackTrackInfo.Converged=1;
    
    
    if J>J0*CtrlVar.NewtonAcceptRatio
        % CtrlVar.BackTrackBeta=1;  % effectivly forces bactracking
        CtrlVar.InfoLevelBackTrack=1000 ; CtrlVar.doplots=1 ;
        Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp
        
        [gamma,J,BackTrackInfo]=BackTracking(slope0,gamma,J0,J,Func,CtrlVar);
        % gamma has changed so I must recalculate J which becomes J0 in the next iteration
        [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
        
    end
    
    GradNorm=norm(dJdp); % grad norm at the end of the iteration step
    
    CtrlVar.NewtonAcceptRatio= max(J/J0,NewtonAcceptRatioMin) ;
    
    % [p,iU,iL]=kk_proj(p,pub,plb);
    % dp(iU)=0 ; dp(iL)=0;

    
    
    if SubOptimality < 0
       fprintf(' ? \n ') 
    end
    
    fprintf("%5i \t %5i \t %7.4g  \t\t %10.5g \t\t  %10.5g \t\t %10.5g \t\t %g \t\t %g \n",Iteration,RunInfo.Inverse.nFuncEval,gamma,J,J/J0,norm(dp)/norm(p+eps),GradNorm,SubOptimality)
    %fprintf("Iteration %i: \t gamma=%-10g \t \t \t J0=%-15.5g \t\t J=%-15.5g \t\t  J1/J0=%-15.5g \t \t |dp|/|p|=%-15.5g \t |dJ/dp|=%-g \n",Iteration,gamma,J0,J,J/J0,norm(dp)/norm(p+eps),GradNorm)
    
    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];
    
    % FindOrCreateFigure('p')  ; PlotMeshScalarVariable(CtrlVar,MUA,p) ; title(sprintf("p at iteration %i",Iteration))
    % FindOrCreateFigure('dp')  ; PlotMeshScalarVariable(CtrlVar,MUA,dp) ; title(sprintf("dp at iteration %i",Iteration))
   
    p=p+gamma*dp;  % do the update
    
    if norm(GradNorm) < eps
        
        fprintf('UaOptimisation: norm of gradient of the objective function smaller than epsilon.\n')
        fprintf('Exciting inverse optimisation step. \n')
        break
    end
    
    if  ~BackTrackInfo.Converged
        fprintf('UaOptimisation: backtrack set in optimisation did not converge.\n')
        fprintf('Exciting inverse optimisation step. \n')
        break
        
     end
    
    
end


return












%%














GradNorm=norm(dJdp);
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





% Determine initial step size:
% If Hessian is defined, use Newton step.
% If CtrlVar.Inverse.InitialLineSearchStepSize is defined use that
%

% if norm(dJdp)<eps
%
%     fprintf('Norm of the gradient of the objective function is less than eps. \n')
%     fprintf('No further inverse iterations needed/possible. \n')
%     return
%
% end

% determine initial search direction and initial step size for line-search.
if ~(isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0)
    gamma=CtrlVar.Inverse.InitialLineSearchStepSize;
    slope0=-dJdp'*dJdp;
else
    if ~(isdiag(Hess) && sum(diag(Hess))==0) && ~strcmpi(CtrlVar.Inverse.DataMisfit.HessianEstimate,'0')
        dp=-Hess\dJdp;
        gamma=1;
        dJdp=-dp;
        slope0=-dJdp'*dJdp;
    else
        slope0=-dJdp'*dJdp;
        gamma1=-0.01*J0/slope0 ; % linear approx
        p1=p-gamma1*dJdp;
        J1=func(p1);
        
        iCount=0 ; % sometimes the initial guess is so small that there is almost no change in the cost function
        % try to increase gamma by factor of 10 until at least 1% change has been generated.
        while (abs(J1-J0)/J1 < 0.01) && iCount<10
            
            gamma1=gamma1*10 ;
            p1=p-gamma1*dJdp;
            J1=func(p1);
            iCount=iCount+1;
        end
        
        gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadradic approx
        if gamma<0 ; gamma=gamma1; end
    end
    
end


dJdpModified=dJdp;

%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \t \t gamma=%-g \n \n',J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,GradNorm,gamma)


Func=@(gamma) func(p-gamma*dJdp);
[J1,~,~,~,~,RunInfo]=Func(gamma);
nFuncEval=1;


if isnan(J1)
    slope0=-dJdp'*dJdp; gamma=-0.01*J0/slope0 ;  % modification on 23 Jan, 2019. Resetting gamma
    J1=Func(gamma);
    nFuncEval=nFuncEval+1;
end

while RunInfo.Forward.uvIterations==0
    
    % the gamma step caused so little change in the model paramters that the previous J0 uv solution was accepted.
    % So increase gamma
    fprintf(" Increasing the stepsize as the previous one caused insufficient changes in model paramters to require a new uv solution.\n")
    fprintf(" gamma increased from %g to %g \n",gamma,gamma*1000)
    gamma=gamma*1000 ;
    [J1,~,~,~,~,RunInfo]=Func(gamma);
    nFuncEval=nFuncEval+1;
    
end

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

fprintf('\n   It   #fEval      J           I          R         |grad|      gamma   #cgUpdate\n')
%fprintf('123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n')

fprintf('%5i  %5i %10g  %10g  %10g  %10g  %10g  %5i\n',It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,GradNorm,gamma,RunInfo.Inverse.ConjGradUpdate)
CtrlVar.GradientUpgradeMethod=CtrlVar.Inverse.GradientUpgradeMethod;

iBackTry=0;

for It=1:CtrlVar.Inverse.Iterations
    
    
    
    [gamma,JgammaNew,BackTrackingInfoVector]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
    nFuncEval=nFuncEval+BackTrackingInfoVector.nFuncEval;
    
    if ~BackTrackingInfoVector.Converged
        
        fprintf(' Line search has stagnated. \n')
        fprintf(' Try resetting step size and using direction of steepest decent. \n')
        slope0=-dJdp'*dJdp; gamma=-0.01*J0/slope0 ;  % modification on 20 Jan, 2019. Resetting gamma
        dJdpModified=dJdp;
        RunInfo.Inverse.ConjGradUpdate=0;
        Func=@(gamma) func(p-gamma*dJdpModified);
        J1=Func(gamma);
        
        
        
        
        [gamma,JgammaNew,BackTrackingInfoVector]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
        nFuncEval=nFuncEval+BackTrackingInfoVector.nFuncEval;
        
        if ~BackTrackingInfoVector.Converged
            fprintf(' Resetting step size to 1 and using steepest decent did not help, now breaking out.\n')
            break
        end
        
    end
    %
    %     if BackTrackingInfoVector.Converged
    %         iBackTry=0;   % After returning from a successful backtracking, update p
    %         p=p-gamma*dJdpModified;
    %         dJdpLast=dJdp;
    %     else
    %         fprintf(' Line search has stagnated,')
    %         iBackTry=iBackTry+1;
    %         if iBackTry==1
    %             fprintf(' try resetting step size to 1.\n')
    %             gamma=1;
    %             continue
    %         else
    %             fprintf(' and resetting step size to 1 did not help, now breaking out.\n')
    %             break
    %         end
    %     end
    
    p=p-gamma*dJdpModified;
    p=kk_proj(p,pub,plb);
    dJdpLast=dJdp;
    % Get new directional derivative
    [J0,dJdp,Hess,fOuts]=func(p);
    nFuncEval=nFuncEval+1; % here J0 and JgammaNew must be (almost) equal
    
    GradNorm=norm(dJdp)/sqrt(numel(dJdp));
    
    % update search direction.
    
    
    fprintf('%5i  %5i %10g  %10g  %10g  %10g  %10g  %5i\n',It+It0,nFuncEval,J0,fOuts.MisfitOuts.I,fOuts.RegOuts.R,GradNorm,gamma,RunInfo.Inverse.ConjGradUpdate)
    
    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];
    
    if CtrlVar.Inverse.StoreSolutionAtEachIteration
        RunInfo.Inverse.p{It+1}=p;
    end
    
    if norm(dJdp) < eps
        norm(dJdp)
        fprintf('UaOptimisation: norm of gradient of the objective function smaller than epsilon.\n')
        fprintf('Exciting inverse optimisation step. \n')
        return
    end
    
    [dJdpModified,RunInfo]=NextGradient(dJdp,dJdpLast,dJdpModified,CtrlVar,RunInfo);
    
    
    
    
    Func=@(gamma) func(p-gamma*dJdpModified);
    % J0=F(0);  This is not needed because I've calculated J0 as I determined the new gradient.
    
    J1=Func(gamma);
    
    if isnan(J1)
        slope0=-dJdp'*dJdp; gamma=-0.01*J0/slope0 ;  % modification on 20 Jan, 2019. Resetting gamma
        J1=Func(gamma);
    end
    
    
end

end
