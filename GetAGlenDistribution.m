function [AGlen,n]=GetAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

[AGlen,n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);

[AGlen,n]=TestAGlenInputValues(CtrlVar,MUA,AGlen,n); 

end