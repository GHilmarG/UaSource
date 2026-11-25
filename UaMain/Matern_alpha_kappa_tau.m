

function [rhoMatern,sigmaMatern,nuMatern]=Matern_alpha_kappa_tau(alphaMatern,kappaMatern,tauMatern)

%% Calculates the Matern parameter triplet:
% 
% rho, sigma, nu 
% 
% from 
% 
% alpha, kappa, and tau 
%
%%


d = 2;                                  % physical dimension
nuMatern=alphaMatern-d/2;            
rhoMatern=sqrt(8*nuMatern)/kappaMatern;      
tau2Matern=tauMatern*tauMatern; 
sigma2Matern = gamma(nuMatern) / (gamma(alphaMatern) * (4*pi)^(d/2) * kappaMatern^(2*nuMatern) * tau2Matern);
sigmaMatern  = sqrt(sigma2Matern);

% (20 Nov 2026) For nu<=0 (alpha=1 in two dimensions, eg the Tikhonov methodology) rho and sigma are not defined: the pointwise
% variance of the continuous field is infinite. They are returned as NaN (previously 0 and Inf).
if nuMatern<=0
    rhoMatern=NaN ; sigmaMatern=NaN ;
end


end