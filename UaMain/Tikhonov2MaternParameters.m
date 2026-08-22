






function [alphaMatern,tauMatern,kappaMatern,sigma2Matern,nuMatern,rhoMatern]=Tikhonov2MaternParameters(ga,gs,Area)

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


alphaMatern=1;
tauMatern=gs/sqrt(2*Area);

if ga==0
    kappaMatern=0;
else
    kappaMatern=ga/(gs+eps(gs));
end

d=2; 
nuMatern=alphaMatern-d/2;

sigma2Matern = gamma(nuMatern) / (gamma(alphaMatern) * (4*pi)^(d/2) * kappaMatern^(2*nuMatern) * tauMatern^2);
rhoMatern=sqrt(8*nuMatern)/kappaMatern;


end