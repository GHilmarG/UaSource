

function [MetricMatrix,QA,QB,QC]=BuildMetricMatrix(CtrlVar,MUA,isA,isB,isC)

nargoutchk(1,4)
narginchk(2,5)

if nargin <3
    [isA,isB,isC] = isABC(CtrlVar);
end

[QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA,isA,isB,isC);

MetricMatrix=blkdiag(QA,QB,QC);


end