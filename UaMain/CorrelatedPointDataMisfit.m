function [J,dJdf,Corr]=CorrelatedPointDataMisfit(CtrlVar,MUA,O,d,dErr,f,CorrPar)
%%
% Misfit between point measurements d and a nodal finite-element field f, with spatially correlated errors (23 November 2026):
%
% $$ J = \frac12 (\mathbf{O}\mathbf{f}-\mathbf{d})^T \mathbf{C}^{-1} (\mathbf{O}\mathbf{f}-\mathbf{d}) , \qquad
%    \mathbf{C} = \mathbf{O}\,\mathbf{Q}_e^{-1}\mathbf{O}^T + \mathbf{D} $$
%
% where Q_e is the precision matrix of a Matern field (alpha=2, standard deviation CorrPar.sigma, range CorrPar.rho) on the mesh, built
% with the SPDE approach as the priors, and D=diag(dErr.^2) the uncorrelated part (nugget). With the Woodbury identity and
% A = O' D^{-1} O, K = Q_e + A (sparse, factorised once and kept),
%
% $$ J = \frac12 \left( \mathbf{r}^T\mathbf{D}^{-1}\mathbf{r} - \mathbf{b}^T\mathbf{K}^{-1}\mathbf{b} \right) , \quad
%    \mathbf{b}=\mathbf{O}^T\mathbf{D}^{-1}\mathbf{r} , \quad
%    \frac{\partial J}{\partial \mathbf{f}} = \mathbf{b} - \mathbf{A}\mathbf{K}^{-1}\mathbf{b} , \quad
%    \frac{\partial^2 J}{\partial \mathbf{f}^2} = \mathbf{A} - \mathbf{A}\mathbf{K}^{-1}\mathbf{A} $$
%
% The Hessian is not formed: Corr returns A and the Cholesky factor of K, and the correction -A K^{-1} A is applied in
% HessianVectorProduct.m (the part A is the Hessian of PointDataMisfit.m with the errors dErr, already in K_qq).
%
%   [J,dJdf,Corr]=CorrelatedPointDataMisfit(CtrlVar,MUA,O,d,dErr,f,CorrPar)
%
% see also: PointDataMisfit.m, Misfit.m, Jq.m, AssembleHessianBlocks.m, HessianVectorProduct.m
%%
persistent Key0 Corr0
w=1./dErr(:).^2;
Key=[MUA.Nnodes size(O,1) nnz(O) full(sum(O*(1:MUA.Nnodes)')) sum(w) sum(w.*(1:numel(w))') CorrPar.sigma CorrPar.rho];
if isempty(Key0) || ~isequal(Key,Key0)
    [kap,tau]=Matern_rho_sigma(2,CorrPar.rho,CorrPar.sigma);
    CV=CtrlVar;
    CV.Inverse.Matern.B.alpha=2;
    CV.Inverse.Matern.B.kappa=kap;
    CV.Inverse.Matern.B.tau=tau;
    CV.Inverse.Regularize.Field="-B-";
    [~,Qe,~]=BuildRegularisationPrecisionMatrices(CV,MUA,false,true,false);
    Qe=(Qe+Qe')/2;
    m=numel(w);
    A=O'*spdiags(w,0,m,m)*O;
    A=(A+A')/2;
    K=Qe+A;
    K=(K+K')/2;
    [L,p]=chol(K,'lower');
    if p~=0
        error("CorrelatedPointDataMisfit:NotPositiveDefinite","Q_e + O'' D^{-1} O is not positive definite.")
    end
    Corr0.A=A;
    Corr0.L=L;
    Corr0.Qe=Qe;
    Key0=Key;
end
Corr=Corr0;
r=O*f(:)-d(:);
wr=w.*r;
b=O'*wr;
z=Corr.L\b;
J=full(r'*wr-z'*z)/2;
if nargout>1
    dJdf=full(b-Corr.A*(Corr.L'\z));
end
end
