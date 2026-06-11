function dIdCbarrier=Calc_dIdCbarrier(CtrlVar,MUA,C)

narginchk(3,3)

if CtrlVar.isBarrierC
    dIdCbarrier=-CtrlVar.muBarrierCmin./(C-CtrlVar.Cmin./2).^2+CtrlVar.muBarrierCmax./(2*CtrlVar.Cmax-C).^2;
else
    dIdCbarrier=C*0;
end

dIdCbarrier=real(dIdCbarrier);
end