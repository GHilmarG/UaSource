function ddIregddC=Calc_dcIregddC(CtrlVar,CC)
    
    [N,M]=size(CC);
    
     if CtrlVar.isRegC
        ddIregddC=CtrlVar.RegCMultiplier*CC;
    else
        ddIregddC=sparse(1:N,1:N);
     end
        
end