function [rho,rhow,g]=TestDensityInputValues(CtrlVar,MUA,rho,rhow,g)



if numel(rho)==1
    
    fprintf(' rho given by user is a scalar. Assuming that rho is same everywhere. \n')
    rho=rho+zeros(MUA.Nnodes,1);
    
end


end