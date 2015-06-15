function dIregdAGlen=Calc_dIregdAGlen(CtrlVar,CAGlen,AGlen,AGlen_prior)
    
    if CtrlVar.isRegAGlen
        dIregdAGlen=CtrlVar.RegAGlenMultiplier*(CAGlen\(AGlen-AGlen_prior));
    else
        dIregdAGlen=AGlen*0;
    end
    dIregdAGlen=real(dIregdAGlen);
end