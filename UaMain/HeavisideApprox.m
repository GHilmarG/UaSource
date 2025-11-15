



function y = HeavisideApprox(k,x,x0)

if nargin==2 ; error('HeavisideApprox: Need at least two arguments \n') ; end

%%
% Smooth approximation to the Heaviside step function using the logistic function, (an example of a sigmoid function).
%
% $$ \frac{dy}{dx}=\frac{1}{1+e^{-2 k (x-x_0)}} $$
%
% The width of the step is about $1/k$, and the limit $k \to \infty$ is (exact) Heaviside step function
%
% $$y \approx 1 \quad \mathrm{if} \quad x > x_0$$
%
% $$y \approx 0 \quad \mathrm{if} \quad x < x_0$$
%
% $$y =1/2 \quad \mathrm{for} \quad x = x_0$$
%
% Note:
%
% The derivative of this function with respect to x is provided by
%
%   DiracDelta(k,x,x0)
%
% and is know as the density of the logistic distribution 
%
% and the integral by 
%
%   SoftPlus(k,x,x0)
% 
% Heaviside:
%
% $$ H=\frac{1}{1+e^{-x}} $$
%
% Dirac delta:
% 
% $$\delta(x)=dH/dx= \frac{1}{\left (e^{-x/2}+e^{x/2} \right )^2} = H(x) (1-H(x)) $$
%
%
% Softplus:
%
% $$ \int H \, dx = \ln (1+e^{x}) $$
%
%
%
% Scaled and offset Heaviside:
%
% $$ H=\frac{1}{1+e^{-2 k (x-x_0) }} $$
%
% Scaled and offset Dirac delta:
% 
% $$\delta(x)=dH/dx= \frac{2k}{(e^{-k (x-x_0) }+e^{k(x-x_0)})^2}  $$
%
%
% Scaled and offset Softplus:
%
% $$ \int H \, dx = \frac{1}{2k} \, \ln (1+e^{2k(x-x_0)}) $$
%
%
% Also
%
%
% $$ H=\frac{1}{1+e^{-2 k (x-x_0) }} = \frac{1}{2} + \frac{1}{2} \tanh( k (x-x_0) )  $$
%
% $$\delta(x)= \frac{2k}{(e^{-k (x-x_0) }+e^{k(x-x_0)})^2}  = \frac{1}{2} \frac{k}{\cosh^2(k(x-x_0))} $$
%
%%

% y=1./(1+exp(-2*k*(x-x0)));  %  The logistic function, an approximation to the Heaviside function 
%  also same as
y=0.5 + 0.5*tanh(k*(x-x0)); % this does not lead to Inf values when calculating the exponent

end