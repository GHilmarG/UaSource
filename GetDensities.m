

function [rho,rhow,g]=GetDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B)


[rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);


if numel(rho)==1
    
    fprintf(' rho given by user is a scalar. Assuming that rho is same everywhere. \n')
    rho=rho+zeros(MUA.Nnodes,1);
    
end

end
