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
%%%



QA=MUA.QA;
QB=MUA.QB;
QC=MUA.QC;

QB=QB+HBobs;

KJpp=blkdiag(QA,QB,QC) ;


end