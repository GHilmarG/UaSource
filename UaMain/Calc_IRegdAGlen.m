function [RegAGlen,dRegdAGlen]=Calc_IRegdAGlen(CtrlVar,MUA,CAGlen,AGlen,AGlen_prior)

narginchk(5,5)

if CtrlVar.isRegAGlen
    
    N=numel(AGlen);
    Ares=(AGlen-AGlen_prior);
    temp=CAGlen\Ares;
    
    RegAGlen=Ares'*temp/(2*N);
    dRegdAGlen=temp/N;
    
    RegAGlen=CtrlVar.RegAGlenMultiplier*RegAGlen;
    dRegdAGlen=CtrlVar.RegAGlenMultiplier*dRegdAGlen;
else
    RegAGlen=0;
    dRegdAGlen=AGlen*0;
end


end
