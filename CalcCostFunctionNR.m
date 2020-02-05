function [UserVar,r,ruv,rl] = CalcCostFunctionNR(UserVar,CtrlVar,MUA,gamma,F,F0,L,l,cuv,dub,dvb,dl)


narginchk(12,12)

if isnan(gamma) ; error(' gamma is nan ') ; end
if ~isreal(gamma) ; error(' gamma is not real ') ; end

F.ub=F.ub+gamma*dub;
F.vb=F.vb+gamma*dvb;
l.ubvb=l.ubvb+gamma*dl;

Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);

%Ruv=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub+gamma*dub,vb+gamma*dvb,uo,vo,AGlen,n,C,m,alpha,rho,rhow,g);


if ~isempty(L)
    frhs=-Ruv-L'*l.ubvb;
    grhs=cuv-L*[F.ub;F.vb];
    
    %frhs=-Ruv-L'*(l+gamma*dl);
    %grhs=cuv-L*[ub+gamma*dub;vb+gamma*dvb];
    
else
    frhs=-Ruv;
    grhs=[];
end

[r,rl,ruv]=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,F0,"-uv-");



end

