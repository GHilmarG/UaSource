function  [p,RunInfo]=QuasiNewtonInversion2(CtrlVar,func,p,RunInfo,MUA)

[J0,dJdp,Hess,fOuts]=func(p);

if numel(RunInfo.Inverse.Iterations)<=1
    RunInfo.Inverse.Iterations(1)=0;
    RunInfo.Inverse.J(1)=J0;
    RunInfo.Inverse.R(1)=fOuts.R;
    RunInfo.Inverse.I(1)=fOuts.I;
    RunInfo.Inverse.StepSize(1)=0;
end

% Determine initial step size:
% If Hessian is defined, use Newton step.
% If CtrlVar.Inverse.InitialLineSearchStepSize is defined use that
%

if norm(dJdp)<eps
   
    fprintf('Norm of the gradient of the objective function is less than eps. \n')
    fprintf('No further inverse iterations needed/possible. \n')
    return
    
end

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
        gamma=-gamma1*slope0/2/((J1-J0)/gamma1-slope0);  % quadradic approx
        if gamma<0 ; gamma=gamma1; end 
    end
    
end


dJdpModified=dJdp;

%%
fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  \t gamma=%-g \n \n',J0,fOuts.I,fOuts.R,gamma)


F=@(gamma) func(p-gamma*dJdp);
J1=F(gamma);


CtrlVar.BacktrackingGammaMin=CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize;
CtrlVar.BackTrackMinXfrac=CtrlVar.Inverse.MinimumRelativelLineSearchStepSize;
CtrlVar.BackTrackMaxIterations=CtrlVar.Inverse.MaximumNumberOfLineSeachSteps;
CtrlVar.InfoLevelBackTrack=CtrlVar.Inverse.InfoLevel;



It0=RunInfo.Inverse.Iterations(end);

fprintf('\n   It          J           I           R       gamma \n')
%fprintf('123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n')
fprintf('%5i\t%10g\t%10g\t%10g\t%10g \n',It0,J0,fOuts.I,fOuts.R,gamma)

CtrlVar.GradientUpgradeMethod=CtrlVar.Inverse.GradientUpgradeMethod;

for I=1:CtrlVar.Inverse.Iterations
    
    
    [gamma,JgammaNew,BackTrackingInfoVector]=BackTracking(slope0,gamma,J0,J1,F,CtrlVar);
    
      
    if ~BackTrackingInfoVector.converged
        fprintf(' Line search has stagnated, breaking out \n')
        break
    end
   
    p=p-gamma*dJdpModified;
    
    dJdpLast=dJdp; 
    [J0,dJdp,Hess,fOuts]=func(p);   % here J0 and JgammaNew must be (almost) equal
    
  
    fprintf('%5i\t%10g\t%10g\t%10g\t%10g \n',I+It0,J0,fOuts.I,fOuts.R,gamma)
    
    
    dJdpModified=NextGradient(dJdp,dJdpLast,dJdpModified,CtrlVar);
        
    % Here I must add modifications to the gradient such as conjugated gradients or
    % BFSG

    
    F=@(gamma) func(p-gamma*dJdpModified);
    J1=F(gamma); % start with previous gamma

    
    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J0];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.I];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];
    
    
    
end

end
