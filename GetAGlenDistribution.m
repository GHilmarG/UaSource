function [UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F,GF)


narginchk(5,5)
nargoutchk(2,2)

N=nargout('DefineAGlenDistribution');


switch N
    
    case 2
        
        [F.AGlen,F.n]=DefineAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,GF);
        
    case 3
        
        [UserVar,F.AGlen,F.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,GF);
        
    otherwise
        
        error('Ua:GetAGlen','DefineAGlenDistribution must return either 2 or 3 output arguments')
        
end

[F.AGlen,F.n]=TestAGlenInputValues(CtrlVar,MUA,F.AGlen,F.n);

end