function [UserVar,r,ruv,rh]=CalcCostFunctionNRuvh(UserVar,CtrlVar,MUA,gamma,du,dv,dh,u,v,h,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,lambda,dlambda)
% CalcCostFunctionNRuvh(CtrlVar,MUA,gamma,du,dv,dh,u,v,h,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,lambda,dlambda);

nargoutchk(4,4)

if nargin~=34
    error(' wrong number of input arguments ')
end

[UserVar,R]=uvhAssembly(UserVar,CtrlVar,MUA,u+gamma*du,v+gamma*dv,h+gamma*dh,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g);


if ~isempty(lambda)
    R=R+L'*(lambda+gamma*dlambda);
end

if nargout==1
    r=ResidualCostFunction(R,F0);
else
    [r,ruv,rh]=ResidualCostFunction(R,F0);
end


end

