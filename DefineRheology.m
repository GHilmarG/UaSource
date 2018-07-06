function [AGlen,n,C,m,rho,rhow,g]=DefineRheology(Experiment,coordinates,connectivity,s,b,h,S,B,Itime,time,CtrlVar)
    
    fprintf('Using default DefineRheology \n')

    Nnodes=size(coordinates,1);
        
    rho=900+zeros(Nnodes,1) ; rhow=1030; g=9.81/1000;
    
    [UserVar,AGlen,n]=DefineAGlenDistribution(Experiment,coordinates,connectivity,s,b,h,S,B,rho,rhow,Itime,time,CtrlVar);
    [UserVar,C,m]=DefineSlipperyDistribution(Experiment,coordinates,connectivity,s,b,h,S,B,rho,rhow,Itime,time,CtrlVar);
    
    
end
