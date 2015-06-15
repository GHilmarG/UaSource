function dIregdC=Calc_dIregdC(CtrlVar,CC,C,C_prior)

if CtrlVar.isRegC
    dIregdC=CtrlVar.RegCMultiplier*(CC\(C-C_prior));
else
    dIregdC=C*0;
end

dIregdC=real(dIregdC);

end