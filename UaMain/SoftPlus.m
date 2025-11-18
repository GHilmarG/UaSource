

function [y,dydx,ddydxx] = SoftPlus(k,x,x0,options)



%%    [y,dydx,ddydxx] = SoftPlus(k,x,x0,options)
%
% Returns a "soft" approximation to the positive part of x.
%
% y=SoftPlus:
%
% $$ y(x) =\mathrm{SoftPlus}(x) = \frac{1}{2k} \, \ln \left ( 1+e^{2k(x-x_0)} \right ) $$
%
% dy/dx =Logistic function, which is a smooth approximation to the Heaviside step-function:
%
% $$ \frac{dy}{dx}=\mathrm{SoftHeaviside}(x) = \frac{1}{1+e^{-2 k (x-x_0)}} $$
%
% ddy/dxx=Smooth approximation to the Dirac delta function: 
%
% $$ \frac{d^2 y}{dx \, dx}=\mathrm{SoftDiracDelta}(x) = \frac{2 k}{\left (e^{k (x-x_0)} +e^{-k (x-x_0)} \right )^2}$$
%
% See HeavisideApprox.m for relationships to the Heaviside step function and the Dirac delta functions.
%
% In the limit of k goes to +infinity the softplus function is 0 for x<=0 and x for x>=0
%
% The derivative of this function is the logistic function, which is a similar type of an approximation to the Heaviside step function
%
% Options:
% 
%   Plot=true     % creates a plot of the outputs, by default Plot=false
%
% Example:
%
%   x=linspace(-2,2,1000) ; k=10 ; x0=0 ;    [y,dydx,ddydxx]=SoftPlus(k,x,x0,Plot=true);
% 
% See also:
%
%   HeavisideApprox.m 
%   DiracDelta.m
%%


arguments
    k     double
    x     double
    x0    double
    options.Plot  logical = false

end

y = log(1+exp(2*k*(x-x0))) /(2*k) ;

I=isinf(y);  % due to overflow errors, for large positive values of x, we have y=inf, instead of y=x
y(I)=x(I);


dydx=[] ;
ddydxx=[];

if nargout > 1

    dydx=HeavisideApprox(k,x,x0) ;

    if nargout ==3
        ddydxx = DiracDelta(k,x,x0);
    end

end


if options.Plot
    
    figure ;
    plot(x,y,DisplayName="SoftPlus") ;
    hold on ;
    if ~isempty(dydx)
    plot(x,dydx,DisplayName="Logistic") ;
    end
    if ~isempty(ddydxx)
    plot(x,ddydxx,DisplayName="Delta") ;
    end
    AL=axis;
    plot([-1/k,1/k],[AL(4)/2 AL(4)/2],"k-",DisplayName="-1/k to 1/k")
    
    xlabel("x") ; ylabel("y") ; 
    title("Softplus, SoftStep, SoftDelta")  ;
    legend();

end

end