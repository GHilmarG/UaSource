function  [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)
%
% func is the function to me minimized
%  p is the paramter set, i.e. func(p)
%
%  Func is func evaluated as a function of stepsize gamma in the direction of
%  the gradient: Func=@(gamma) func(p-gamma*dJdp);
%

% 
%
narginchk(8,8)
nargoutchk(3,3)

%%
p=p(:);
p=kk_proj(p,pub,plb);

[J,dJdp,Hess,fOuts]=func(p); RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;



GradNorm=norm(dJdp)/sqrt(numel(dJdp));
RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
RunInfo.Inverse.J=[RunInfo.Inverse.J;J];
RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;NaN];




fprintf("   It \t #fEval\t   gamma   \t\t       J \t\t\t\t  J/J0 \t\t\t     |dp|/|p| \t\t |dJ/dp| \t sub-obtimality gap \t theta (degrees) \n")
fprintf("%5i \t %5i \t %7.4g  \t\t %10.5g \t\t  %10.5g \t\t %10.5g \t\t %10.5g \t %10.5g \t %10.5g\n",0,RunInfo.Inverse.nFuncEval,NaN,J,NaN,NaN,GradNorm,NaN,NaN)



NewtonAcceptRatioMin=0.8 ;
CtrlVar.NewtonAcceptRatio=NewtonAcceptRatioMin ;
CtrlVar.BackTrackMinXfrac=1e-3 ;
gamma=NaN;




for Iteration=1:CtrlVar.Inverse.Iterations
    
    J0=J;
    
%     Zeros=sparse(MUA.Nnodes,MUA.Nnodes);
%     pH=[M Zeros ; ...
%         Zeros M ] ;
%     lambda=1;
%     lambda=lambda*mean(abs(diag(Hess)))/mean(abs(diag(pH)));
%     dp=-(Hess+lambda*pH)\dJdp ;   % Newton system
%     
    dp=-Hess\dJdp ;   % Newton system

    slope0=dJdp'*dp;
    SubOptimality=-dJdp'*dp/2  ;  % sub-optimality at the beginning of the iteration step

    % A sensible step length is something like
    if isnan(gamma)
        gamma=-0.01*J0/slope0  ;
        % gamma=-dJdp'*dp/(dp'*Hess*dp);  % or if I beleve in the Hessian then...
        
        % gamma*dp/norm<1e-4  -> gamma=1e-4 norm(p)/norm(dp)
        % gamma=1e-4*norm(p)/norm(dp);
        
    end
    
    theta=real(acos(-dJdp'*dp/(norm(dJdp)*norm(dp))));
    % [~,theta2]=FE_inner_product(-dJdp,dp,MUA.M) ; % This is actually the correct angle, but as far as I can see the difference is small
    % [real(theta) real(theta2)]*180/pi
    
    if CtrlVar.InfoLevelInverse>=10
        
        
        fprintf("angle between search direction and steepest-decent direction is %f degrees.\n",theta*180/pi)
        
        if CtrlVar.InfoLevelInverse>=1000
            %% Test quadratic model and slope
            
            CtrlVar.InfoLevelBackTrack=1000 ; % also get info on backtracking
            
            Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp
            
            f=[] ; x=[];
            [f,x]=TestQuadradicModel(f,x,Func,Hess,dp,gamma,J0,dJdp,slope0);
            TestSlope(Func,0,0.01*gamma,slope0) ;
            
        end
        
        % JQuad=J0+gamma*dJdp'*dp+gamma^2*dp'*Hess*dp/2 ; % must do this before I calculate the new Hessian
    end
    
 
    
    [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;

    while (isnan(J) || any(isnan(dJdp))) && ~(any(isnan(p))  || any(isnan(dp)))  
         warning('UaOptimisationHessianBased:NaNinObjectiveFunction',' Objective function returned a NaN. Trying a new step size.') ;
         gamma=gamma/100 ; [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
    end
    
    BackTrackInfo.Converged=1;
    if J>J0*CtrlVar.NewtonAcceptRatio
        
        Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp
        
        [gamma,J,BackTrackInfo]=BackTracking(slope0,gamma,J0,J,Func,CtrlVar);
        % gamma has changed so I must recalculate J which becomes J0 in the next iteration
        [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
        
    end
    
    GradNorm=norm(dJdp)/sqrt(numel(dJdp));
    CtrlVar.NewtonAcceptRatio= max(J/J0,NewtonAcceptRatioMin) ;

    fprintf("%5i \t %5i \t %7.4g  \t\t %10.5g \t\t  %10.5g \t\t %10.5g \t\t %10.5g \t %10.5g \t\t %10.5g \n",...
        Iteration,RunInfo.Inverse.nFuncEval,gamma,J,J/J0,norm(dp)/norm(p+eps),GradNorm,SubOptimality,theta*180/pi)
    
    
    RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+1];
    RunInfo.Inverse.J=[RunInfo.Inverse.J;J];
    RunInfo.Inverse.R=[RunInfo.Inverse.R;fOuts.RegOuts.R];
    RunInfo.Inverse.I=[RunInfo.Inverse.I;fOuts.MisfitOuts.I];
    RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNorm];
    RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gamma];
    
    
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


end

