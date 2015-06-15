
function IRegC=Calc_IRegC(CtrlVar,CC,C,C_prior)
    
    
    if CtrlVar.isRegC
        Cres=(C-C_prior);
        IRegC=CtrlVar.RegCMultiplier*Cres'*(CC\Cres)/2    ;
    else
        IRegC=0;
    end
    
end
