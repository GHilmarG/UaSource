function dIregdAGlen=Calc_dIregdAGlen(CtrlVar,MUA,CAGlen,AGlen,AGlen_prior)

narginchk(5,5)

if CtrlVar.isRegAGlen
    dIregdAGlen=CtrlVar.RegAGlenMultiplier*(CAGlen\(AGlen-AGlen_prior));
else
    dIregdAGlen=AGlen*0;
end
dIregdAGlen=real(dIregdAGlen);
end