

function [G,QA,QB,QC]=BuildMetricMatrix(CtrlVar,MUA,Meas)

nargoutchk(1,4)
narginchk(3,3)

[isA,isB,isC] = isABC(CtrlVar);

[QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA,isA,isB,isC);

if isB

    if ~isempty(Meas.Bobs)
   
        O=Meas.BO;
        nMeas=numel(Meas.Bobs);
        iSigma=sparse(1:nMeas,1:nMeas,1./Meas.BErr.^2,nMeas,nMeas);
        HBobs=O' * iSigma * O;
        QB=QB+HBobs;
    end

end



G=blkdiag(QA,QB,QC);


end