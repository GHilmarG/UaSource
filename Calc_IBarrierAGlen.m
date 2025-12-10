function IBarrierAGlen=Calc_IBarrierAGlen(CtrlVar,AGlen)
    
    if CtrlVar.isBarrierAGlen
        
        %IBarrierMin=-CtrlVar.muBarrierAGlenmin*sum(log(AGlen-CtrlVar.AGlenmin/2));
        %IBarrierMax=-CtrlVar.muBarrierAGlenmax*sum(log(2*CtrlVar.AGlenmax-AGlen));
        
        IBarrierMin=CtrlVar.muBarrierAGlenmin*sum(1./(AGlen-CtrlVar.AGlenmin/2));
        IBarrierMax=CtrlVar.muBarrierAGlenmax*sum(1./(2*CtrlVar.AGlenmax-AGlen));
        
        
        IBarrierAGlen=IBarrierMin+IBarrierMax;
        if isinf(IBarrierAGlen) ;
            warning('MisfitFunction:NaN',' IBarrierAGlen is inf')
            IBarrierAGlen=0;
        end
    else
        IBarrierAGlen=0;
    end
    
end
