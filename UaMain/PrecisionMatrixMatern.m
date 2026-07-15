



function Q=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,ga,gs,Methodology)

narginchk(7,7)
nargoutchk(1,1)


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

Area=MUA.Area;
M=MUA.M;
D=MUA.Dxx+MUA.Dyy;

% I discovered that it can happen that the smallest eigenvalue of D is slightly
% negative!!! This must be due to numerical rounding errors when assembling Dxx and Dyy. I for example found a case where the
% two smallest eigenvalues of Dyy were -1.14405445408737e-16 and  -8.99887803162969e-17. One approach of dealing with this
% would be to add eps to the diagonal of Dxx and Dyy. This should not really be an issue, unless there is next-to-no
% M contribution being added.

Ieps=sparse(1:MUA.Nnodes,1:MUA.Nnodes,eps);
D=D+Ieps ; 

%  MUA.M=MassMatrix2D1dof(MUA);
% [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);


if Methodology=="-Tikhonov-" 

    %Q=(gsB.^2.*(Dxx+Dyy)+gaB.^2.*M)/Area;

    alphaMatern=1;
    tauMatern=gs/sqrt(Area);
    if ga==0
        kappaMatern=0;
    else
        kappaMatern=ga/(gs+eps(gs));
    end

end

A = kappaMatern^2*M + D;

if alphaMatern==1

    % Q1
    Q = tauMatern^2 * A ;

elseif alphaMatern==2

    % Q2
    row_sums = sum(M,2);
    Q = tauMatern^2 * A * sparse(1:n,1:n,1./row_sums,n,n) * A;

else

    error("alpha must be either equal to 1 or 2.\n")

end




end

