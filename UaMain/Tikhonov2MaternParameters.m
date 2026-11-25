






function [alphaMatern,tauMatern,kappaMatern]=Tikhonov2MaternParameters(ga,gs,Area)

nargoutchk(3,3)
narginchk(3,3)


%% Calculates the Matern parameters from the Tikhonov parameters
%
%
% Matern formulation is:
%
% $$Q_{\alpha=1} = \tau^2\left(\kappa^2 M + D\right) $$
%
% using Tikhonov formulation I've written the precision matrix as:
%
% $$Q = \frac{1}{\mathcal{A}} ( \gamma_a^2 M + \gamma_s^2 D ) $$
%
% (20 Nov 2026) The regularisation term is R = 1/2 dp' Q dp (Regularisation.m), and the Tikhonov term is defined as
%
% $$R = \frac{1}{2\mathcal{A}} \left( \gamma_a^2 \, dp^T M dp + \gamma_s^2 \, dp^T D dp \right) $$
%
% so that Q = (gamma_a^2 M + gamma_s^2 D)/Area. Previously Q was taken to be (gamma_a^2 M + gamma_s^2 D)/(2 Area), ie with tau a factor
% sqrt(2) too small, so that the regularisation was half as strong as defined.
%
% therefore
%
% $$\alpha=1 $$
%
% $$\tau=\frac{\gamma_s}{\sqrt{\mathcal{A} }} $$
%
% $$\kappa=\frac{\gamma_a}{\gamma_s} $$
%
%%


if isapprox(ga,0)
    fprintf("The amplitude regularisation term in the Tikhonov formulation is close to zero. \n")
    fprintf("Run continues, but precision matrices likley to be singular or near singular.  \n")
    warning("Tikhonov2MaternParameters:gaTooSmall","The amplitude regularisation term in the Tikhonov formulation is close to zero. \n")
end

alphaMatern=1;
tauMatern=gs/sqrt(Area);   % (20 Nov 2026) was gs/sqrt(2*Area), see above

if ga==0
    kappaMatern=0;
else
    kappaMatern=ga/(gs+eps(gs));
end



end