function IBarrierC=Calc_IBarrierC(CtrlVar,C)
    
    if CtrlVar.isBarrierC
        
        %IBarrierMin=-CtrlVar.muBarrierCmin*sum(log(C-CtrlVar.Cmin/2));
        %IBarrierMax=-CtrlVar.muBarrierCmax*sum(log(2*CtrlVar.Cmax-C));
        
        
        IBarrierMin=CtrlVar.muBarrierCmin*sum(1./(C-CtrlVar.Cmin/2));
        IBarrierMax=CtrlVar.muBarrierCmax*sum(1./(2*CtrlVar.Cmax-C));
        
        
        IBarrierC=IBarrierMin+IBarrierMax;
        if isinf(IBarrierC) ;
            warning('MisfitFunction:inf',' IBarrierC is inf')
            IBarrierC=0;
        end
    else
        IBarrierC=0;
    end
    
end

