
function KFpp=Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)




if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    KFAA=FAA_v2(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    error("Fpp:BnotImplemented","Fpp not yet implemented for B")
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    KFCC=FCC_v2(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
end


%% The forward model only has terms involving either A or C, never both. So there are no cross-terms the the Fpp matrix is a diagonal block matrix
% Note that this is not the case when expanding the inversion to B. Then we have terms involving A, C and h, where h=s-B.

if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true) && contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    KFpp = blkdiag(KFAA,KFCC);
elseif contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    KFpp=KFAA;
elseif contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    KFpp=KFCC;
else
    error("PsiTimesddFuvdpdp:FellThrough","Fell through")
end


end

