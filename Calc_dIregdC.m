function dIregdC=Calc_dIregdC(CtrlVar,MUA,CC,C,C_prior)

narginchk(5,5)

if CtrlVar.isRegC
    dIregdC=CtrlVar.RegCMultiplier*(CC\(C-C_prior));
else
    dIregdC=C*0;
end

dIregdC=real(dIregdC);

end