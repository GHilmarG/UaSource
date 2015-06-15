function [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,ubvbLambda,udvdLambda,dIdu,kv,rh,nlInfo]=...
    CalcMisfitFunction(Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,AGlen,C,n,m,alpha,rho,rhow,g,GF,Priors,Meas)

narginchk(26,26)

[ub,vb,ud,vd,ubvbLambda,udvdLambda,kv,rh,nlInfo]= uv(CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,AGlen,C,n,m,alpha,rho,rhow,g,GF);
[J,Idata,IRegC,IRegAGlen,dIdu,IBarrierC,IBarrierAGlen]=MisfitFunction(Experiment,CtrlVar,MUA,ub,vb,ud,vd,AGlen,C,Priors,Meas);

end

