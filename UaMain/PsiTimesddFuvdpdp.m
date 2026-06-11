

function HPsiddFdpdp=PsiTimesddFuvdpdp(CtrlVar,MUA,F,uAdjoint,vAdjoint)


if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    HPsiddFdAdA=PsiTimesddFuvdAdA(CtrlVar,MUA,F,uAdjoint,vAdjoint);
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    error("PsiTimesddFuvdpdp:BnotImplemented","PsiTimesddFuvdBdB not yet implemented")
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    HPsiddFdCdC=PsiTimesddFuvdCdC(CtrlVar,MUA,F,uAdjoint,vAdjoint);
end

if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true) && contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    HPsiddFdpdp = blkdiag(HPsiddFdAdA,HPsiddFdCdC);
elseif contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    HPsiddFdpdp=HPsiddFdAdA;
elseif contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    HPsiddFdpdp=HPsiddFdCdC;
else
    error("PsiTimesddFuvdpdp:FellThrough","Fell through")
end


end

