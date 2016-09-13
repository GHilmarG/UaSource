function [UserVar,r,ruv,rh,rl]=CalcCostFunctionNRuvh(UserVar,CtrlVar,MUA,gamma,dub,dvb,dh,ub,vb,h,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,l,dl,cuvh)


nargoutchk(5,5)
narginchk(35,35)


[UserVar,R]=uvhAssembly(UserVar,CtrlVar,MUA,ub+gamma*dub,vb+gamma*dvb,h+gamma*dh,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g);


if ~isempty(L)
    frhs=-R-L'*(l+gamma*dl);
    grhs=cuvh-L*[ub+gamma*dub;vb+gamma*dvb;h+gamma*dh];
else
    frhs=-R;
    grhs=[];
end

[r,rl,ruv,rh]=ResidualCostFunction(frhs,grhs,F0,MUA.Nnodes);



end

