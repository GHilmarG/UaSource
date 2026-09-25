function KJpp=Jpp(CtrlVar,MUA)

narginchk(2,2)

%% Builds the Hessian term, Jpp, which is
%
% $$\delta^2_{pp} J $$ 
%
% This is the Hessian resulting from explicit dependency of J on p. 
%
% This is therefore a very easy term to calculate.
%
%% Contributions
%
% There are two:
%
% # the regularisation terms, whose Hessians are the precision matrices QA, QB and QC
% # the misfit with respect to direct observations of B, if any
%
% Both are built in
%
%   [G,QA,QB,QC,HBobs]=BuildMetricMatrix(CtrlVar,MUA,Meas)
%
% and stored as fields of MUA. Note that BuildMetricMatrix.m deliberately keeps HBobs out of QB, and out of the
% metric G. That is the correct split: the metric used for the Riesz mapping and for the trust region is the prior
% precision alone, whereas the Hessian needs the prior plus the data term. The two are therefore added together here
% and not in BuildMetricMatrix.m
%
% The scaling matches without any rescaling. BobsMisfit.m in Regularisation.m forms
%
% $$ J_{B_{obs}} = \frac{1}{2} \, r^T \Sigma^{-1} r , \qquad r = O B - B_{obs} $$
%
% so that its Hessian is $O^T \Sigma^{-1} O$, which is exactly HBobs.
%
% Inactive fields drop out automatically: BuildRegularisationPrecisionMatrices.m returns empty matrices for the
% fields not being inverted for, and blkdiag then omits those blocks. HBobs is likewise empty unless B is inverted
% for and direct observations of B have been provided.
%
%  see also: BuildMetricMatrix.m, BuildRegularisationPrecisionMatrices.m, Fpp.m, CalcDirectAdjointHessian.m
%
%%%



QA=MUA.QA;
QB=MUA.QB;
QC=MUA.QC;

% add the Hessian of the direct-B-observation misfit to the B block
%
% The isfield test is needed because MUA.HBobs is only set where BuildMetricMatrix.m has been called, i.e. within
% InvertForModelParameters.m

if isfield(MUA,"HBobs") && ~isempty(MUA.HBobs)
    QB=QB+MUA.HBobs;
end


KJpp=blkdiag(QA,QB,QC) ;


end