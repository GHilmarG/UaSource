function KJpp=Jpp(CtrlVar,MUA,Meas)

[isA,isB,isC] = isABC(CtrlVar) ; 

QA=MUA.QA;
QB=MUA.QB;
QC=MUA.QC;

if isB

    % This does not change during the inversion, so I could do this once.
    % However, this is fast, and furthermore I'm expecting to have to change this in the near future
    Berr=full(sqrt(spdiags(Meas.BCov)));
    iBerr = spdiags(1./Berr,0,MUA.Nnodes,MUA.Nnodes);
    ddRdBmeasBmeas = iBerr*MUA.M*iBerr/MUA.Area;

    QB=QB+ddRdBmeasBmeas;


end


KJpp=blkdiag(QA,QB,QC) ;



end