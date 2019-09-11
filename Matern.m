function [r,nu,kappa,sigma2]=Matern(sigma,alpha,rho,dimention,distance)

%%
%
%   sigma    :   marginal variance
%   alpha    :   alpha/2 is the exponent in the fractional Helmholtz eqaution 
%    rho     :   distance where correlation falls to 0.1
%  dimention : spatial dimention
%
%
% alpha=nu+d/2
% nu=alpha-d/2
%
% for d=2 (two spatial dimentions) a
% nu=2-1=1
%
%  rho=sqrt(8 nu)/kappa
%  then given rho
%  kappa=sqrt(8 nu)/rho
%
% Example :
%
%   d=2; alpha=2 ; rho=1e3 ; sigma=10 ; nu=NaN ; kappa=NaN ; dist=linspace(1,2000,100) ;
%   [r,nu,kappa]=Matern(sigma,alpha,rho,d,dist); 
%   figure ; plot(dist,r)
%
%%


% this gives a correlation of about 0.1 at the distance r
nu=alpha-dimention/2;  % ie nu=1 for alpha=2 and dimention=2
kappa=sqrt(8*nu)/rho;

sigma2=gamma(nu) /  ( gamma(nu+dimention/2)*(4*pi)^(dimention/2)*kappa^(2*nu));

x=kappa * distance ;

r = sigma2 * x.^nu .* besselk(nu,x)  / (2^(nu-1)*gamma(nu)) ;

end


