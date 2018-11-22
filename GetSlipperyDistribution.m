function [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F)

narginchk(4,4)
nargoutchk(2,2)


N=nargout('DefineSlipperyDistribution');


switch N
    
    case 2
        
        [F.C,F.m]=DefineSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
        
    case 3
        
        [UserVar,F.C,F.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
        
    otherwise
        
        error('Ua:GetSlipperyDistribution','DefineSlipperyDistribution must return either 2 or 3 arguments')
        
end



[F.C,F.m]=TestSlipperinessInputValues(CtrlVar,MUA,F.C,F.m);
[F.C,iU,iL]=kk_proj(F.C,CtrlVar.Cmax,CtrlVar.Cmin);


end




