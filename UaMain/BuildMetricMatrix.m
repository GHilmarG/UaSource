

function [G,QA,QB,QC,HBobs]=BuildMetricMatrix(CtrlVar,MUA,Meas)

nargoutchk(1,5)
narginchk(3,3)

HBobs=[]; 
[isA,isB,isC] = isABC(CtrlVar);

[QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA,isA,isB,isC);

if isB

    if ~isempty(Meas.Bobs)

        O=Meas.BO(Meas.BInside,:);
        Inside=Meas.BInside;
        Bobs=Meas.Bobs(Inside);
        BErr=Meas.BErr(Inside);
        nMeas=numel(Bobs);

        iSigma=sparse(1:nMeas,1:nMeas,1./BErr.^2,nMeas,nMeas);
        HBobs=O' * iSigma * O;
        HBobs=0.5*(HBobs+HBobs');
        % QB=QB+HBobs;
    end
end



G=blkdiag(QA,QB,QC);


end