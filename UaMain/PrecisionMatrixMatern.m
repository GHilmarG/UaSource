



function PrecisionMatrixMatern(CtrlVar,M,D,alpha,kappa,tau)

%%
% *Precision matrix*
%
% $$ Q_{\alpha=1} = \tau^2\left(\kappa^2 M + K\right) $$
%
% $$Q_{\alpha=2} = \tau^2 \; (\kappa^2 M + K)\; M^{-1} \; (\kappa^2 M + K)$$
%
% But the inverse of the mass matrix is, in general, dense. So we simply cannot use this approach without a modification.
%
% Instead we do:
%
% $$Q_{\alpha=2} = \tau^2 \; (\kappa^2 M + K)\; \tilde{M}^{-1} \; (\kappa^2 M + K)$$
%
% where $$\tilde{M}$$ is a diagonal row-sum lumping
%
% $$\tilde{M}_{ik} = \delta_{ik} \, \sum_j M_{ij} $$
%
% According to Lindgre-Rue-Lindström, this does not lead to a new leading-term error in the FE solution, and is absorbed into
% the existing discretization error budget.
%
%
%
% *Mat\'ern Covariance*
%
% If we solve the equation:
%
% $$\tau \; (\kappa^2 - \Delta)^{\alpha/2} \;  B(x) = W(x) $$
%
% $$x \in \mathcal{R}^d $$
%
% with $$W(x)$$ uncorrelated Gaussian white noise. Then the resulting covariance is
%
% $$ C(r) = \frac{\sigma^2}{2^{\nu - 1}\, \Gamma(\nu)} \;  \left(\kappa \|r\|\right)^{\nu} \; K_{\nu}\!\left(\kappa
% \|r\|\right) $$
%
% where
%
% $$r=(x-x_0,y-y_0) $$ as the covariance is isotropic. consequences
%
%
% Here we have
%
% $$\nu = \alpha - \frac{d}{2}$$
%
% where $$d$$ is the dimension, i.e $$d=2$$ in the xy plane, and the marginal variance is
%
% $$\sigma^2 = \frac{\Gamma(\nu)}{\Gamma(\alpha)\, (4\pi)^{d/2}\, \kappa^{2\nu}\, \tau^2} $$
%
%
% *Priors*
%
% the expression for the prior, is on the form
%
%
% $$-\log P(B) = \frac{1}{2}(B-B_{prior})^{T} Q (B-B_{prior}) + \frac{1}{2}\, \log \, \left|Q^{-1}\right| + \mathrm{const} $$
%
% see also: Matern.m which gives an more detailed information about the Matern process, and provides the link between the
% Matern hyper-parameters and the $$\kappa$$ and $$\tau$$ parameters.
%%

%QA=0.5*(gsA.^2.*(Dxx+Dyy)+gaA.^2.*M)/Area; % This is the precision matrix
NB=(gsB.^2.*(Dxx+Dyy)+gaB.^2.*M)/Area;
RB=dpB'*NB*dpB/2;               %       R: Regularisation term for B (a scalar)
dRdB=(NB*dpB).*dBfactor;        %   dR/dB:  (a vector)
ddRdBB=NB.*dBfactor;            % exact, or simply the correct, Hessian of the regularization term
% To do: I could add "RHB=E" to CtrlVar.Inverse.Hessian. Right now I do the exact (E) Hessian evaluation here.

if alpha==1



elseif alpha==2


else

    error("alpha must be either equal to 1 or 2.\n")

end




end

