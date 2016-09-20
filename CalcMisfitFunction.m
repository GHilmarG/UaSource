function [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
    CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,Priors,Meas)



narginchk(25,25)

            
[UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,ubvbL]=...
    uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);


[J,Idata,IRegC,IRegAGlen,dIdu,IBarrierC,IBarrierAGlen]=MisfitFunction(UserVar,CtrlVar,MUA,ub,vb,ud,vd,AGlen,C,Priors,Meas);

end

