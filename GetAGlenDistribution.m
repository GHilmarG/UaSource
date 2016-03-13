function [AGlen,n]=GetAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

[AGlen,n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);


if numel(AGlen)==1
    
    fprintf(' AGlen given by user is a scalar. Assuming that AGlen is same everywhere. \n')
    if  CtrlVar.AGlenisElementBased
        AGlen=AGlen+zeros(MUA.Nele,1);
    else
        AGlen=AGlen+zeros(MUA.Nnodes,1);
    end
    
end

end