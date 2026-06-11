function dRegCdC=Calc_dIregdC(CtrlVar,MUA,CC,C,C_prior)

% In the future, should get rid of this and use [RegC,dRegCdC]=Calc_IRegC(CtrlVar,CC,C,C_prior)
% instead which calculates both the value and the gradient

narginchk(5,5)

if CtrlVar.isRegC
    N=numel(C);
    Cres=(C-C_prior);
    dRegCdC=CtrlVar.RegCMultiplier*(CC\Cres)/N;
else
    dRegCdC=C*0;
end

dRegCdC=real(dRegCdC);

end