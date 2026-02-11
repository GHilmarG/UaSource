

function [y,dydx,ddydxx] = SoftPlus(k,x,x0,options)



%%    [y,dydx,ddydxx] = SoftPlus(k,x,x0,options)
%
% Returns a "soft" approximation to the positive part of x-x0.
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

if numel(x0)>1
    y(I)=x(I)-x0(I);
else
    y(I)=x(I)-x0;
end

dydx=[] ;
ddydxx=[];

if nargout > 1

    dydx=HeavisideApprox(k,x,x0) ;

    %dydx=1./(1+exp(-2*k*(x-x0)));

    if nargout ==3
        SoftDiracDelta = DiracDelta(k,x,x0);
        ddydxx=SoftDiracDelta;
        % ddydxx=sparse(1:numel(SoftDiracDelta),1:numel(SoftDiracDelta),SoftDiracDelta,numel(SoftDiracDelta),numel(SoftDiracDelta));
    end



end


if options.Plot


    xmx0=x-x0;
    [xSorted,I]=sort(xmx0);
    %ddydxx=diag(ddydxx);

    FindOrCreateFigure("SoftPlus")

    T=tiledlayout("flow");

    T1=nexttile;

    plot(xmx0(I),y(I),DisplayName="Soft-Plus",LineWidth=2,Color="k") ;
    xlim([-10/k 10/k])
    xlabel("$x-x_0$",Interpreter="latex") ;
    ylabel("$y$",Interpreter="latex") ;
    title("Soft-Plus")
    lg1=legend(Location="best");

    hold on ;
    if ~isempty(dydx)
        T2=nexttile;
        plot(xmx0(I),dydx(I),LineWidth=2,DisplayName="Soft-Heaviside",color="r") ;
        xlim([-10/k 10/k])
        xlabel("$x-x_0$",Interpreter="latex") ;
        title("Soft Heaviside")
        ylabel("$y$",Interpreter="latex") ;
        lg2=legend(Location="best");

    end
    if ~isempty(ddydxx)
        T3=nexttile;
        plot(xmx0(I),ddydxx(I),DisplayName="Soft-Delta",LineWidth=2,Color="b") ;
        xlim([-10/k 10/k])
        xlabel("$x-x_0$",Interpreter="latex") ;
        ylabel("$y$",Interpreter="latex") ;
        lg2=legend(Location="best");
        title("Soft delta")

    end
    % AL=axis;
    % plot([-1/k,1/k],[AL(4)/2 AL(4)/2],"k-",DisplayName="-1/k to 1/k")

    xlabel("$x-x_0$",Interpreter="latex") ;


    lg=legend(Location="best");

end

end