



function [R2,x,lambda,dx,dlambda,Slope0,gammamin,BackTrackInfo,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,c,H0,R20,K0,R0,g0,h0,KK0)


exitflag=0 ; 


[dx,dlambda]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);
Slope0=2*R0'*K0*dx ;


if Slope0 > 0
    fprintf("lsqStepUa: Exiting because slope at origin in line search positive (Slope=%g) \n",Slope0)
    R2=R20 ; g2=g20 ; x=x0 ; lambda=lambda0; dx=0 ; dlambda=0 ; gammamin=nan ; BackTrackInfo=[]; 
    exitflag=1 ; 
    return
end



gammaEst=-R0'*K0*dx/(dx'*(KK0)*dx) ;

funcBackTrack=@(gamma) R2func(gamma,dx,dlambda,fun,x0,lambda0) ;


R2=nan;

CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;

%CtrlVar.NewtonAcceptRatio=0.001; 


[gammamin,R2min,BackTrackInfo]=BackTracking(Slope0,gammaEst,R20,R2,funcBackTrack,CtrlVar);

dx=gammamin*dx  ;
dlambda=gammamin*dlambda ;

x=x0+dx ;
lambda=lambda0+dlambda ;

R2=R2min ;


end