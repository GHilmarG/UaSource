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
% Note: This does not yet have the B misfit term related to direct observations of B 
%
% This is already done in 
%
%   [G,QA,QB,QC,HBobs]=BuildMetricMatrix(CtrlVar,MUA,Meas)
%
% so should be easy to add, and should already be a field of MUA
%
%%%



QA=MUA.QA;
QB=MUA.QB;
QC=MUA.QC;


KJpp=blkdiag(QA,QB,QC) ;


end