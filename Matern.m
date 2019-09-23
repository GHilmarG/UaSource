function [r,nu,kappa,sigma2,Cov,Realisation]=Matern(alpha,rho,dimention,distance,sigma2,DoPlots)

%%
%
%   alpha    :   alpha/2 is the exponent in the fractional Helmholtz eqaution 
%    rho     :   distance where correlation falls to 0.1
%  dimention : spatial dimention
%
%
% alpha=nu+d/2
% nu=alpha-d/2
% sigma2          :   marginal variance
%
% for d=2 (i.e. two spatial dimentions) we have nu=2-1=1
%
%  rho=sqrt(8 nu)/kappa
%  then given rho
%  kappa=sqrt(8 nu)/rho
%
%
% If sigma2 is not given on input, it is calculated based on expression in:
%   Lindgren, F., Rue, H., & Lindström, J. (2011). 
%
% Note that if sigma2 is specified, the variance of the Helmoltz equation with
% the specified alpha, rho and dimention will still be sigma2 as given by Eq 1
% below!  Only input sigma2 if you interested in returning and using the
% covariance matrix Cov and the Realisation.
%
%
% Example :
%
%   d=2; alpha=2 ; rho=1e3 ; sigma=10 ; nu=NaN ; kappa=NaN ; dist=linspace(1,2000,100) ;
%   [r,nu,kappa]=Matern(sigma,alpha,rho,d,dist); 
%   figure ; plot(dist,r)
%
%%

if nargin<6
    DoPlots=false;
end

Cov=[] ; Realisation=[] ; 

% this gives a correlation of about 0.1 at the distance r
nu=alpha-dimention/2;  % ie nu=1 for alpha=2 and dimention=2

kappa=sqrt(8*nu)/rho;

if nargin<5
    % Eq 1
    sigma2=gamma(nu) /  ( gamma(nu+dimention/2)*(4*pi)^(dimention/2)*kappa^(2*nu));
end

x=kappa * distance ;

% this sigma2 should be referred to as sigma^2
% some other sources appear to use a different defintion of rho 
% for example rho in https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
% appear to be 1/2 of the rho I use here
% My notation is based on: 
%   Lindgren, F., Rue, H., & Lindström, J. (2011). 
%   An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), 423–498. https://doi.org/10.1111/j.1467-9868.2011.00777.x
r = real(sigma2 * x.^nu .* besselk(nu,x)  / (2^(nu-1)*gamma(nu)) );


%%

if nargout>4 || DoPlots
    F = griddedInterpolant(distance,r) ;
    D=ndgrid(distance,distance) ;
    D=abs(D-D') ;
    
    Cov=F(D) ;
    R=sqrtm(Cov);
    Realisation=R*randn(numel(distance),1) ;
    
    if DoPlots
        figure;
        for I=1:3
            y=R*randn(numel(distance),1) ;
            plot(distance/1000,y) ; ylabel('Matern realisation') ; xlabel('distance (km)')
            hold on
            fprintf(' Expected variance %f  \t  estimated variance %f \n ',sigma2,var(y))
        end
        title(sprintf(' A few examples of Matern realisations with rho=%i and sigma=%i',rho,sigma2))
    end
    
end

end






