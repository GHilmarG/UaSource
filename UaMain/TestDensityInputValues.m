function [rho,rhow,g]=TestDensityInputValues(CtrlVar,MUA,rho,rhow,g)



if isscalar(rho)
    
    fprintf(' rho given by user is a scalar. Assuming that rho is spatially constant. \n')
    rho=rho+zeros(MUA.Nnodes,1);
    
end


errorStruct.identifier = 'Ua:NaNinInput';

if any(isnan(rho))
    errorStruct.message = 'nan in rho';
    error(errorStruct)
end


if any(isnan(rhow))
    errorStruct.message = 'nan in rhow';
    error(errorStruct)
end




end