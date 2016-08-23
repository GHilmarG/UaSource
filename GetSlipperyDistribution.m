function [UserVar,C,m]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)



nOut=nargout;

if nOut~=3
    error('Ua:GetSlipperyDistribution','Need 3 output arguments')
end

N=nargout('DefineSlipperyDistribution');


switch N
    
    case 2
        
        [C,m]=DefineSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
        
    case 3
        
        [UserVar,C,m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
        
    otherwise
        
        error('Ua:GetSlipperyDistribution','DefineSlipperyDistribution must return either 2 or 3 arguments')
        
end



[C,m]=TestSlipperinessInputValues(CtrlVar,MUA,C,m);


end




