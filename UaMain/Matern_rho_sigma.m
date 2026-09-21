function [kappaMatern,tauMatern,nuMatern]=Matern_rho_sigma(alphaMatern,rhoMatern,sigmaMatern)

%% Calculates the Matern parameters kappa, tau, nu from alpha, rho and sigma
%
% Exact inverse of Matern_alpha_kappa_tau.m
%
%   rho   : correlation range (distance at which correlation ~0.1)
%   sigma : marginal standard deviation of the field
%   alpha : smoothness index, must be > d/2 = 1
%
% $$\nu=\alpha-d/2, \quad \kappa=\sqrt{8\nu}/\rho, \quad
%   \tau=\frac{1}{\sigma}\sqrt{\frac{\Gamma(\nu)}{\Gamma(\alpha)(4\pi)^{d/2}}}\,\kappa^{-\nu}$$
%

narginchk(3,3)

d = 2;                                   % physical dimension
nuMatern = alphaMatern - d/2;

if nuMatern <= 0
    error("Matern_rho_sigma:alpha", ...
        "alpha must exceed d/2=%g for a valid Matern field, but alpha=%g.",d/2,alphaMatern)
end

kappaMatern = sqrt(8*nuMatern)/rhoMatern;

tau2Matern  = gamma(nuMatern) / ...
    (gamma(alphaMatern)*(4*pi)^(d/2)*kappaMatern^(2*nuMatern)*sigmaMatern^2);
tauMatern   = sqrt(tau2Matern);

end
