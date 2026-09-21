

function [QA,QB,QC]=BuildRegularisationPrecisionMatrices(CtrlVar,MUA,isA,isB,isC)

nargoutchk(3,3)
narginchk(5,5)

%% Builds a metric matrix
%
% The metric matrix $G$ plays the role of mimicking the effect of an inner product between two functions $f$ and $g$ such
% that
%
% $$ \langle f(x) \vert g(x) \rangle = \mathbf{x}^{t} G \mathbf{y}$$
% 
%
% Here I'm using the exact same matrices as I would in the regularization for logA, B, and logC. These are the precision
% matrices used in the cost terms when regularizing logA, B and logC.
%
% This seems to give the right limit behavior. With
%
% $$J=M(p)+R(p) $$
%
% and
%
% $$\nabla^2 R= Q $$
%
% setting $G=Q$ gives the Newton step
%
% $$ G^{-1} \nabla J = Q^{-1} g_{M} + (p-p_{\mathrm{prior}}) $$
%
% i.e. the $G$-steepest-descent direction takes exact Newton step on the regularization term and prior-preconditions the
% misfit gradient.
%
% This is dimensional consistent. The logA, B, and logC blocks have different physical dimensions. But the QA, QB, QC
% matrices have the the corresponding inverse units because each carries its own $\tau^2$ in the Matern covariance. (This
% only holds if one uses the Matern covariance.) 
%
% However there is no $A/C$ coupling and this is a weakness. Not sure how to deal with this though...
%
%
% Both $\alpha=1$ and $\alpha=3$ will work, but the former gives a degenerate random field, and the latter arguably a too
% smooth one.
%
% 
%
%%




QA=[] ; QB=[]; QC=[];

if isA
    alphaMatern=CtrlVar.Inverse.Matern.logAGlen.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logAGlen.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logAGlen.tau;
    gsA=CtrlVar.Inverse.Regularize.logAGlen.gs;
    gaA=CtrlVar.Inverse.Regularize.logAGlen.ga;
    QA=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaA,gsA,CtrlVar.Inverse.Methodology);
end


if isB
    alphaMatern=CtrlVar.Inverse.Matern.B.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.B.kappa;
    tauMatern=CtrlVar.Inverse.Matern.B.tau;
    gsB=CtrlVar.Inverse.Regularize.B.gs;
    gaB=CtrlVar.Inverse.Regularize.B.ga;
    QB=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaB,gsB,CtrlVar.Inverse.Methodology);
end


if isC
    alphaMatern=CtrlVar.Inverse.Matern.logC.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logC.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logC.tau;
    gsC=CtrlVar.Inverse.Regularize.logC.gs;
    gaC=CtrlVar.Inverse.Regularize.logC.ga;

    QC=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaC,gsC,CtrlVar.Inverse.Methodology);
end


end

