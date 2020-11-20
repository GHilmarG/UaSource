function [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F)
    
    narginchk(4,4)
    nargoutchk(2,2)
    
    
    N=nargout('DefineSlipperyDistribution');
    
    
    if any(CtrlVar.SlidingLaw==["Budd","W-N0"]) && N<4
        error("GetSlipperyDistribution:nargout","When using Budd sliding law, DefineSlipperyDistribution.m must return 4 parameters [UserVar,C,m,q] ")
    end
    
    
    if any(CtrlVar.SlidingLaw==["Tsai","Coulomb","Cornford","Umbi","W","W-N0","minCW-N0","C","rpCW-N0","rCW-N0"])  && N<5
        fprintf("\n \n When using the sliding law %s, DefineSlipperyDistribution.m must return 5 parameters [UserVar,C,m,q,muk]. \n",CtrlVar.SlidingLaw)
        fprintf(" Note: The sliding law only depends on the parameters C, m, muk. (so you can set, for example, q=NaN.)  \n")
        error("GetSlipperyDistribution:nargout","Incorrect number of output parameters in DefineSlipperiness ")
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




