






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
% $$Q = \frac{1}{2 \mathcal{A}} ( \gamma_a^2 M + \gamma_s^2 D ) $$
%
% therefore
%
% $$\alpha=1 $$
%
% $$\tau=\frac{\gamma_s}{\sqrt{2 \mathcal{A} }} $$
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
tauMatern=gs/sqrt(2*Area);

if ga==0
    kappaMatern=0;
else
    kappaMatern=ga/(gs+eps(gs));
end



end