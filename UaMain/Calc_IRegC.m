
function [RegC,dRegCdC]=Calc_IRegC(CtrlVar,MUA,CC,C,C_prior)

narginchk(5,5)

if CtrlVar.isRegC
    
    % R= (C-C_prior)' CC^{-1} (C-C_prior)  / (2N)
    
    N=numel(C);
    Cres=(C-C_prior);
    temp=CC\Cres;
    
    RegC=Cres'*temp/(2*N)   ;
    dRegCdC=temp/N;
    
    RegC=CtrlVar.RegCMultiplier*RegC ;
    dRegCdC=CtrlVar.RegCMultiplier*dRegCdC;
    
else
    RegC=0;
    dRegCdC=C*0;
end

end
