




function Hessian=HessianABC(p,lambda,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)

narginchk(12,12)

CtrlVar.JGH.CalcHessian=true; % But will only do so if the number of output arguments is also 3 or greater
[~,~,Hessian]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ; 





end