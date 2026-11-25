

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
        % QB=QB+HBobs;   % HBobs is deliberately kept out of QB, as QB is also used to evaluate R and dR/dB in Regularisation.m
    end
end



% The metric is the Hessian of the regularisation term, and optionally also includes the Hessian of the misfit with respect to the
% direct measurements of B, HBobs (in the metric only, see above). (14 November 2026)
if ~isempty(HBobs) && isfield(CtrlVar.Inverse,"Metric") && isfield(CtrlVar.Inverse.Metric,"IncludeBobsHessian") ...
        && CtrlVar.Inverse.Metric.IncludeBobsHessian
    G=blkdiag(QA,QB+HBobs,QC);
else
    G=blkdiag(QA,QB,QC);
end


end