function [r,ruv,rl] = CalcCostFunctionNR...
    (gamma,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,MUA,alpha,rho,rhow,g,F0,L,l,dl,CtrlVar,cuv)


narginchk(24,24)

if isnan(gamma) ; error(' gamma is nan ') ; end
if ~isreal(gamma) ; error(' gamma is not real ') ; end

if CtrlVar.CalvingFrontFullyFloating
    R=RTF(s,S,B,h,ub+gamma*dub,vb+gamma*dvb,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar);
else
    R=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub+gamma*dub,vb+gamma*dvb,AGlen,n,C,m,alpha,rho,rhow,g);
end

R2=0;
if ~isempty(l)
    
    
    R=R+L'*(l+gamma*dl);
    R2=cuv-L*[ub+gamma*dub;vb+gamma*dvb];  % (2016 Sept)
    
end



%fgamma=ResidualCostFunction(R,F0);

ruv=R'*R;
rl=R2'*R2;
ruv=real(full(ruv)) ; rl=real(full(rl));

%fprintf(' r1=%g \t r2=%g \n',r1,r2)

r=full((ruv+rl)/(F0'*F0));
r=real(r);




end

