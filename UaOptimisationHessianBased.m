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
TrustRadius=NaN ; 
for Iteration=1:CtrlVar.Inverse.Iterations
    
    J0=J;
    
    
    dp=-Hess\dJdp ;   % Newton system
    
    slope0=dJdp'*dp;
    
    % A sensible step length is something like
    if isnan(gamma)
        gamma=-0.01*J0/slope0  ;
        % gamma=-dJdp'*dp/(dp'*Hess*dp);  % or if I beleve in the Hessian then...
    end
    
    %% Test quadratic model and slope 
    
    
    Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp

    f=[] ; x=[]; 
    [f,x]=TestQuadradicModel(f,x,Func,Hess,dp,gamma,J0,dJdp,slope0);
    TestSlope(Func,0,0.01*gamma,slope0) ;
    %%
    
    if isnan(TrustRadius)
        TrustRadius=0.01*norm(dp) ;
    end
    
    SubOptimality=-dJdp'*dp/2  ;  % sub-optimality at the beginning of the iteration step
    
    
    
    theta=acos(-dJdp'*dp/(norm(dJdp)*norm(dp)));
    fprintf("angle between search direction and steepest-decent direction is %f degrees.\n",theta*180/pi)
    
    JQuad=J0+gamma*dJdp'*dp+gamma^2*dp'*Hess*dp/2 ; % must do this before I calculate the new Hessian
    [J,dJdp,Hess,fOuts]=func(p+gamma*dp);  RunInfo.Inverse.nFuncEval=RunInfo.Inverse.nFuncEval+1;
    
    
    fprintf('J=%g \t JQuad=%g \t (J-JQuad)/J=%g \n',J,JQuad,(JQuad-J)/J)
    
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

  
    
    fprintf("%5i \t %5i \t %7.4g  \t\t %10.5g \t\t  %10.5g \t\t %10.5g \t\t %g \t\t %g \n",Iteration,RunInfo.Inverse.nFuncEval,gamma,J,J/J0,norm(dp)/norm(p+eps),GradNorm,SubOptimality)
    
    
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
    
    
    % The quadratic approximation is:
    % j=J(0)+dJdp'*dp+dp'*H*dp/2 ;
    % This is not the Cauchy point because I'm not enforcing norm(gamma dp) = gamma norm(dp) < TrustRadius
    rhok=(J-J0)/(gamma*dJdp'*dp+gamma^2*dp'*Hess*dp/2) ;
    fprintf(' Ratio between actual and predicted reduction of the cost function based on the quadratic model is: %g \n',rhok)
    
    etav=0.9 ; etas=0.1 ; gi=2 ; gd=0.5 ; 
    if rhok>etav  % very succesfull
        gammaOld=gamma ; 
        gammaNew=-dJdp'*dp/(dp'*Hess*dp); 
        gamma=min(gammaNew,gi*gammaOld) ; 
        TrustRadius=TrustRadius*gi;
        fprintf(' Trust radius: very successful.\n')
        fprintf(' Selecting new gamma based on the Cauchy point. gamma= %g \n',gamma)
        
        
    elseif rhok>etas
        fprintf(' Trust radius: successful.\n')
    else
        fprintf(' Trust radius: not successful.\n')
        TrustRadius=TrustRadius*gd; 
    end
    fprintf(' Trust radius=%g \n',TrustRadius)
    
    
end


return


end

