

function [MetricMatrix,QA,QB,QC]=BuildMetricMatrix(CtrlVar,MUA,isA,isB,isC)

nargoutchk(1,4)
narginchk(5,5)


[QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA,isA,isB,isC);

MetricMatrix=blkdiag(QA,QB,QC);


end