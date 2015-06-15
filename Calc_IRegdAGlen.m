function IRegAGlen=Calc_IRegdAGlen(CtrlVar,CAGlen,AGlen,AGlen_prior)
    
    if CtrlVar.isRegAGlen
        Ares=(AGlen-AGlen_prior);
        IRegAGlen=CtrlVar.RegAGlenMultiplier*Ares'*(CAGlen\Ares)/2;
    else
        IRegAGlen=0;
    end
    
    
end
