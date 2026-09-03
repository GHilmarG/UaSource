

function MetricMatrix=BuildMetricMatrix(CtrlVar,MUA)

 [QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA);



MetricMatrix=blkdiag(QA,QB,QC);


end