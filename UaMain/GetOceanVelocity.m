function [UserVar,uo,vo]=GetOceanVelocity(UserVar,CtrlVar,MUA,BCs,ub,vb,ud,vd,uo,vo,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m)

[UserVar,uo,vo]=DefineOceanVelocity(UserVar,CtrlVar,MUA,CtrlVar.time,ub,vb,ud,vd,uo,vo,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
                                   
end