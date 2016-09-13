function [UserVar,r,ruv,rl] = CalcCostFunctionNR...
    (UserVar,CtrlVar,MUA,gamma,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,l,dl,cuv)


narginchk(25,25)

if isnan(gamma) ; error(' gamma is nan ') ; end
if ~isreal(gamma) ; error(' gamma is not real ') ; end


R=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub+gamma*dub,vb+gamma*dvb,AGlen,n,C,m,alpha,rho,rhow,g);


if ~isempty(L)
    frhs=-R-L'*(l+gamma*dl);
    grhs=cuv-L*[ub+gamma*dub;vb+gamma*dvb];
else
    frhs=-R;
    grhs=[];
end

[r,rl,ruv]=ResidualCostFunction(frhs,grhs,F0,MUA.Nnodes);



end

