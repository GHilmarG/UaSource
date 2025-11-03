function [rho,rhow,g]=TestDensityInputValues(CtrlVar,MUA,rho,rhow,g)



if numel(rho)==1
    
    fprintf(' rho given by user is a scalar. Assuming that rho is same everywhere. \n')
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