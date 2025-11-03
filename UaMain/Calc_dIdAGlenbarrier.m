function dIdAGlenbarrier=Calc_dIdAGlenbarrier(CtrlVar,MUA,AGlen)

narginchk(3,3)

if CtrlVar.isBarrierAGlen
    
    %dIdAGlenbarrier=-CtrlVar.muBarrierAGlenmin./(AGlen-CtrlVar.AGlenmin/2)+CtrlVar.muBarrierAGlenmax./(2*CtrlVar.AGlenmax-AGlen); % log
    
    dIdAGlenbarrier=-CtrlVar.muBarrierAGlenmin./(AGlen-CtrlVar.AGlenmin/2).^2+...
        CtrlVar.muBarrierAGlenmax./(2*CtrlVar.AGlenmax-AGlen).^2;
    
else
    dIdAGlenbarrier=AGlen*0;
end

dIdAGlenbarrier=real(dIdAGlenbarrier);
end
