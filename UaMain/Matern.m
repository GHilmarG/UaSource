





function [r,nu,kappa,sigma2Helmholtz,Cov,Realisation]=Matern(alpha,rho,d,x,sigma2,DoPlots)

%%
%
% Remark: Noticed in July 2026, that this was done for $$\tau=1$$, which is OK, but not the most general expression
%
% If we solve the equation:
%
% $$\tau \; (\kappa^2 - \Delta)^{\alpha/2} \;  B(x) = W(x) $$
%
% $$x \in \mathcal{R}^d $$ 
%
% with $$W(x)$$ uncorrelated Gaussian white noise. Then the resulting covariance is 
%
% $$ C(h) = \frac{\sigma^2}{2^{\nu - 1}\, \Gamma(\nu)} \left(\kappa \|h\|\right)^{\nu} K_{\nu}\!\left(\kappa  \|h\|\right) $$
%
% where
%
% $$h=(x-x_0,y-y_0) $$ as the covariance is isotropic. 
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
% The correlation function $$\rho(h)$$ is simply the covariance normalized to give $\rho(h=0)=1$$ and is
%
% $$ \rho(h) = \frac{C(h)}{\sigma^2} $$
% 
%
% We have several parameters that need to be selected. Of those $$d$$ is not really a free parameter but given by
% the physical dimension of the problem. 
%
% I can select $$\kappa$$, $$\alpha$$ and $$\tau$$ with $$d=2$$ and then set
%
%
% $$\nu = \alpha - \frac{d}{2}$$
% 
% and 
%
% $$\sigma^2 = \frac{\Gamma(\nu)}{\Gamma(\alpha)\, (4\pi)^{d/2}\, \kappa^{2\nu}\, \tau^2} $$
%
% This makes sense if I'm thinking about the parameters of the differential operator as independent. 
%
% Alternatively, I can select $$\sigma^2$$ and $\alpha$ with $$d=2$$, but this is not enough and I will additionally
% need to fix $$\kappa$$. This is often done by selecting the distance over which the correlation drops to 0.1 which
% turns out to be related to $$\kappa$$ and $\nu$$, so this is really just another way of selecting $$\kappa$$. This
% gives
%
% $$\kappa=\frac{\sqrt{8 \nu}}{\rho}$$
%
% where $$\rho$$ in this expression is the distance at which the correlation drops to 0.1
% 
% This is arguably the most physical way of thinking about this.I select $$\alpha$$ which controls the smoothness. This
% gives me $$\nu$$ since $$d$ is not really a free parameter. Then fix a correlation length $$\rho$$ which gives me
% $$\kappa$$ from the above equation. I now also need to select $$\sigma^2$$, which relates to the uncertainty in the
% magnitude of the prior. I can then calculate $$\tau$$ from the expression for $$\sigma^2$$.
%
%
% 
%  [r,nu,kappa,sigma2,Cov,Realisation]=Matern(alpha,rho,d,x,sigma2,DoPlots)
%
% Calculates the Matern covariance defined as
%
% $$ r(x)=\frac{\sigma^2}{2^{\nu-1} \Gamma(\nu)} (\kappa \, x )^{\nu} K_{\nu}(\kappa x)$$
%
% Inputs:
%
%    alpha   :   alpha/2 is the exponent in the fractional Helmholtz equation
%    rho     :   distance where correlation falls to 0.1
%    d       :   spatial dimension
%
%    sigma2  :   marginal variance, (optional input)
%   DoPlots  :   plots of realizations if set to true (optional input)
%
% Outputs: 
%
%       nu              :   nu=alpha-d/2; 
%    kappa              :   the wave-number in the Helmholtz equation
%    sigma2Helmholtz    :   marginal variance of the Helmholtz equation for the given
%                           kappa and alpha values 
%    Cov                :   covariance matrix with the marginal covariance
%                           sigma2, if sigma2 is provided, otherwise with the
%                           marginal covariance sigma2Helmholtz
%    Realisation        :   One realization of a Matern process. 
% 
% $$\nu=\alpha-d/2$$
%
% for d=2 (i.e. two spatial dimensions) we have nu=2-1=1
%
%  rho=sqrt(8 nu)/kappa then given rho 
%
% $$\kappa=\frac{\sqrt{8 \nu}}{\rho}$$
%
%
% If sigma2 is not given on input, it is calculated based on expression in:
%   Lindgren, F., Rue, H., & Lindström, J. (2011).
%
% Note that if $\sigma^2$ is specified, the marginal variance of the Helmoltz equation with
% the specified alpha, rho and dimention will still be sigma2 as given 
%
%
% $$ \sigma^2=\frac{\Gamma(\nu) }{\Gamma(\nu+d/2) (4 \pi)^{d/2} \kappa^{2 \nu} } $$ 
%
% 
% by the equation above!
% 
% Only input sigma2 if you interested in returning and using the
% covariance matrix Cov and the Realisation.
%
%
% Example :
%
%   d=2; alpha=2 ; rho=4e3 ; 
%   dist=linspace(0,1e4,1e3) ; [r,nu,kappa,sigma2]=Matern(alpha,rho,d,dist);
%   figure ; plot(dist,r) ; xlabel('distance') ; ylabel('Matern')
%   title(['$$ \sigma^2= $$',num2str(sigma2),' $$ \rho=$$',num2str(rho)],'interpreter','latex')
%%

if nargin<6
    DoPlots=false;
end

Cov=[] ; Realisation=[] ; 

% this gives a correlation of about 0.1 at the distance r
nu=alpha-d/2;  % ie nu=1 for alpha=2 and dimension=2

kappa=sqrt(8*nu)/rho;

% sigma^2 : the marginal variance 
sigma2Helmholtz=gamma(nu) /  ( gamma(alpha)*(4*pi)^(d/2)*kappa^(2*nu));

if nargin<5 || isempty(sigma2)
    % Eq 1
    sigma2=sigma2Helmholtz ; 
end

dist=x ; % for plotting purposes
x=kappa * x ;

% this sigma2 should be referred to as sigma^2
% some other sources appear to use a different definition of rho 
% for example rho in https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
% appear to be 1/2 of the rho I use here
% My notation is based on: 
%   Lindgren, F., Rue, H., & Lindström, J. (2011). 
%   An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), 423–498. https://doi.org/10.1111/j.1467-9868.2011.00777.x

r = real(sigma2 * x.^nu .* besselk(nu,x)  / (2^(nu-1)*gamma(nu)) );

% the expression above for r fails for x=0, although the limit x^nu beselk(nu,x)  for x  to 0 is well defined, in fact
%
%   lim x to 0 of x^nu besselk(nu,x) = 1 
%

Ix0=x==0 ; 
r(Ix0) = real(sigma2  / (2^(nu-1)*gamma(nu)) );

%%

if nargout>4 || DoPlots
    F = griddedInterpolant(x,r) ;
    D=ndgrid(x,x) ;
    D=abs(D-D') ;
    
    Cov=F(D) ;
    R=sqrtm(Cov);
    Realisation=R*randn(numel(x),1) ;
    
    if DoPlots
        figure;
        for I=1:3
            y=R*randn(numel(x),1) ;
            plot(dist/1000,y) ; ylabel('Matern realisation') ; xlabel('distance (km)')
            hold on
            fprintf(' Expected variance %f  \t  estimated variance %f \n ',sigma2,var(y))
        end
        title(sprintf(' A few examples of Matern realisations with rho=%i and sigma=%i',rho,sigma2))
    end
    
end

end






