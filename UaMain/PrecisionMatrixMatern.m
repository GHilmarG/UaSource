



function PrecisionMatrix(CtrlVar,M,D)

%%
%
% $$\tau \; (\kappa^2 - \Delta)^{\alpha/2} \;  B(x) = W(x) $$
%
% $$x \in \mathcal{R}^d $$ 
%
% with $$W(x)$$ uncorrelated Gaussian white noise.
%
% $$ C(h) = \frac{\sigma^2}{2^{\nu - 1}\, \Gamma(\nu)} \left(\kappa \|h\|\right)^{\nu} K_{\nu}\!\left(\kappa  \|h\|\right) $$
%
% where
%
% $$h=(x-x_0,y-y_0) $$ as the covariance is isotropic. 
%
% $$\nu = \alpha - \frac{d}{2}$$
%
% $$\sigma^2 = \frac{\Gamma(\nu)}{\Gamma(\alpha)\, (4\pi)^{d/2}\, \kappa^{2\nu}\, \tau^2} $$
%
% $$ (\kappa^2 - \Delta)^{\alpha/2}\big(\tau\, B(x)\big) = \mathcal{W}(x)$$
%
%%

%QA=0.5*(gsA.^2.*(Dxx+Dyy)+gaA.^2.*M)/Area; % This is the precision matrix
NB=(gsB.^2.*(Dxx+Dyy)+gaB.^2.*M)/Area;
RB=dpB'*NB*dpB/2;               %       R: Regularisation term for B (a scalar)
dRdB=(NB*dpB).*dBfactor;        %   dR/dB:  (a vector)
ddRdBB=NB.*dBfactor;            % exact, or simply the correct, Hessian of the regularization term
% To do: I could add "RHB=E" to CtrlVar.Inverse.Hessian. Right now I do the exact (E) Hessian evaluation here.

end

