



function [R2,x,lambda,dx,dlambda]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,c,H0,R20,K0,R0,g0,h0,KK0)




nx=numel(x0);



[dx,dlambda]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);


Slope0=2*R0'*K0*dx ;
gammaEst=-R0'*K0*dx/(dx'*(KK0)*dx) ;

funcBackTrack=@(gamma) R2func(gamma,dx,dlambda,fun,x0,lambda0) ;

CtrlVar.BacktrackIteration=nan  ;
R2=nan;

CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ; CtrlVar.NewtonAcceptRatio=0.001; CtrlVar.doplots=1 ;
[gammamin,R2min,BackTrackInfo]=BackTracking(Slope0,gammaEst,R20,R2,funcBackTrack,CtrlVar);

dx=gammaminNewton*dx  ;
dlambda=gammamin*dlambdaN ;

x=x0+dx ;
lambda=lambda0+dlambda ;

R2=R2min ; 


end