function [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F)
    
    narginchk(4,4)
    nargoutchk(2,2)
    
    
    N=nargout('DefineSlipperyDistribution');
    
    
    if CtrlVar.SlidingLaw=="Budd" && N<4
        error("GetSlipperyDistribution:nargout","When using Budd sliding law, DefineSlipperyDistribution.m must return 4 parameters [UserVar,C,m,q] ")
    end
    
    
    if CtrlVar.SlidingLaw=="Tsai" && N<5
        error("GetSlipperyDistribution:nargout","When using Tsai sliding law, DefineSlipperyDistribution.m must return 5 parameters [UserVar,C,m,q,muk] ")
    end
    
    
    switch N
        
        case 2
            
            [F.C,F.m]=DefineSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
        case 3
            
            [UserVar,F.C,F.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
        case 4
            
            [UserVar,F.C,F.m,F.q]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
        case 5
            
            [UserVar,F.C,F.m,F.q,F.muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
            
        otherwise
            
            error('Ua:GetSlipperyDistribution','DefineSlipperyDistribution must return either 2 or 3 arguments')
            
    end
    
    
    
    [F.C,F.m,F.q,F.muk]=TestSlipperinessInputValues(CtrlVar,MUA,F.C,F.m,F.q,F.muk);
    [F.C,iU,iL]=kk_proj(F.C,CtrlVar.Cmax,CtrlVar.Cmin);
    
    
end




