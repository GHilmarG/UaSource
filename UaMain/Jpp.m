function KJpp=Jpp(CtrlVar,MUA)

narginchk(2,2)

QA=MUA.QA;
QB=MUA.QB;
QC=MUA.QC;



KJpp=blkdiag(QA,QB,QC) ;



end