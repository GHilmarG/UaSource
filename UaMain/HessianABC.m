




function Hessian=HessianABC(p,lambda,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)


%
% fmincon with the interior-point method requires the Hessian being returned in a separate function. 
%
%   hessian = hessianfcn(x,lambda)
%
% where lambda is a structure with the Lagrange multiplier associated with the nonlinear constraints.
%
% I don't use any non-linear constraints in the A/, B, C inversions, so lambda is not used, but still required as an input to
% this function.
%
%
%

narginchk(12,12)

CtrlVar.JGH.CalcHessian=true; % But will only do so if the number of output arguments is also 3 or greater
[~,~,Hessian]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ; 





end