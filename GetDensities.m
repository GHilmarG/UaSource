

function [rho,rhow,g]=GetDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B)


[rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);

[rho,rhow,g]=TestDensityInputValues(CtrlVar,MUA,rho,rhow,g);



end
