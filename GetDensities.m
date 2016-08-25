

function [UserVar,rho,rhow,g]=GetDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)


nOut=nargout;
if nOut~=4
    error('Ua:GetDensities','Need 4 output arguments')
end

N=nargout('DefineDensities');

switch N
    
    case 3
        
        [rho,rhow,g]=DefineDensities(CtrlVar.Experiment,CtrlVar,MUA,time,s,b,h,S,B);
        
    case 4
        
        [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B);
        
    otherwise
        
        error('Ua:GetDensities','Need 3 or 4 outputs')
        
end

[rho,rhow,g]=TestDensityInputValues(CtrlVar,MUA,rho,rhow,g);

end
