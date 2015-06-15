function fgamma = CalcCostFunctionNR(gamma,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,MUA,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar)
                 

   if isnan(gamma) ; error(' gamma is nan ') ; end
   if ~isreal(gamma) ; error(' gamma is not real ') ; end
    
    if CtrlVar.CalvingFrontFullyFloating
        R=RTF(s,S,B,h,ub+gamma*dub,vb+gamma*dvb,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar);
    else
        R=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub+gamma*dub,vb+gamma*dvb,AGlen,n,C,m,alpha,rho,rhow,g);
    end
    
    if ~isempty(lambda)
        R=R+L'*(lambda+gamma*dlambda);
    end
    
   
    
    fgamma=ResidualCostFunction(R,F0);
   
    
    
end

