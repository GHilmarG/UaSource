function y = HeavisideApprox(k,x,x0)

if nargin==2 ; error('HeavisideApprox: Need three arguments \n') ; end

%%
% Approximation to the Heaviside step  function
%  The width of the step is about $1/k$
%  the limit $k \to \infty$ is Heaviside step function
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
%
%%


y=1./(1+exp(-2*k*(x-x0)));  %  The logistic function
%  also same as
%  0.5 + 0.5 tanh(x/2)

end