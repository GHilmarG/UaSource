

function Q=HessPosDefAddition(CtrlVar,MUA)


QA=[] ; QB=[] ; QC=[];

% A
if contains(lower(CtrlVar.Inverse.Regularize.Field),'logaglen')

    % regularize log10(AGlen)

    gsA=CtrlVar.Inverse.Regularize.logAGlen.gs;
    gaA=CtrlVar.Inverse.Regularize.logAGlen.ga;
    alphaMatern=CtrlVar.Inverse.Matern.logAGlen.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logAGlen.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logAGlen.tau;


    QA=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaA,gsA,CtrlVar.Inverse.Methodology);

end


% C
if contains(lower(CtrlVar.Inverse.InvertFor),'logc')  % this includes both c and logc inversion

    gsC=CtrlVar.Inverse.Regularize.logC.gs;
    gaC=CtrlVar.Inverse.Regularize.logC.ga;
    alphaMatern=CtrlVar.Inverse.Matern.logC.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logC.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logC.tau;

    QC=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaC,gsC,CtrlVar.Inverse.Methodology);


end



% B
if contains(CtrlVar.Inverse.InvertFor,'-B-')

    gsB=CtrlVar.Inverse.Regularize.B.gs;
    gaB=CtrlVar.Inverse.Regularize.B.ga;
    alphaMatern=CtrlVar.Inverse.Matern.B.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.B.kappa;
    tauMatern=CtrlVar.Inverse.Matern.B.tau;
    QB=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaB,gsB,CtrlVar.Inverse.Methodology);


end


Q=blkdiag(QA,QB,QC);


end

