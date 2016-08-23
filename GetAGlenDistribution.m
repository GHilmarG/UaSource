function [UserVar,AGlen,n]=GetAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)



nOut=nargout;
if nOut~=3
    error('Ua:AGlen','Need 3 output arguments')
end

N=nargout('DefineAGlenDistribution');


switch N
    
    case 2
        
        [AGlen,n]=DefineAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
        
    case 3
        
        [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
        
    otherwise
        
        error('Ua:GetAGlen','DefineAGlenDistribution must return either 2 or 3 output arguments')
        
end

[AGlen,n]=TestAGlenInputValues(CtrlVar,MUA,AGlen,n);

end