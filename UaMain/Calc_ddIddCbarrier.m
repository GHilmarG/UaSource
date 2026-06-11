function ddIddCbarrier=Calc_ddIddCbarrier(CtrlVar,MUA,C)

narginchk(3,3)

N=length(C);

if CtrlVar.isBarrierC
    
    ddIddCbarrier=2*CtrlVar.muBarrierCmin./(C-CtrlVar.Cmin./2).^3+2*CtrlVar.muBarrierCmax./(2*CtrlVar.Cmax-C).^3;
    
else
    
    ddIddCbarrier=C*0;
    
end

ddIddCbarrier=sparse(1:N,1:N,ddIddCbarrier);

end