function dIregdAGlen=Calc_dIregdAGlen(CtrlVar,MUA,CAGlen,AGlen,AGlen_prior)

% In the future, get rid of this and use
%[RegAGlen,dRegdAGlen]=Calc_IRegdAGlen(CtrlVar,CAGlen,AGlen,AGlen_prior)
% to calclate both value and gradient.


narginchk(5,5)

if CtrlVar.isRegAGlen
    N=numel(AGlen);
    Ares=(AGlen-AGlen_prior);
    dIregdAGlen=(1/N)*CAGlen\Ares;
    dIregdAGlen=CtrlVar.RegAGlenMultiplier*dIregdAGlen;
else
    dIregdAGlen=AGlen*0;
end
dIregdAGlen=real(dIregdAGlen);

end
