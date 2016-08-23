function [UserVar,ub,vb,ud,vd]=GetStartVelValues(UserVar,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m)


nOut=nargout;
if nOut~=5
    error('Ua:GetStartVelValues','Need 5 output arguments')
end

N=nargout('DefineStartVelValues');

switch N
    
    case 4
        
        
        [ub,vb,ud,vd]=DefineStartVelValues(CtrlVar.Experiment,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
        
        
    case 5
        
        [UserVar,ub,vb,ud,vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
        
end

end

