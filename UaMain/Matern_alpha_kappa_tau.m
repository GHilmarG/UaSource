

function [rhoMatern,sigmaMatern,nuMatern]=Matern_alpha_kappa_tau(alphaMatern,kappaMatern,tauMatern)

%% Calculates the Matern parameters rho, sigma, nu from alpha, kappa, and tau 



d = 2;                                  % physical dimension
nuMatern=alphaMatern-d/2;            
rhoMatern=sqrt(8*nuMatern)/kappaMatern;      
tau2Matern=tauMatern*tauMatern; 
sigma2Matern = gamma(nuMatern) / (gamma(alphaMatern) * (4*pi)^(d/2) * kappaMatern^(2*nuMatern) * tau2Matern);
sigmaMatern  = sqrt(sigma2Matern);


end