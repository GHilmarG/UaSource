function [f,x]=TestQuadradicModel(fIn,xIn,Func,Hess,dp,gamma,J0,dJdp,slope0,xMax,N)

if nargin<10 || isempty(xMax)
    gammaMin=-dJdp'*dp/(dp'*Hess*dp);  % this should actually always be equal to 1, if dp is the Newton step
    xMax=min(gammaMin,2*gamma) ; 
end

if nargin<11 || isempty(N)
    N=8;
end


Nquad=200 ; xQuad=linspace(0,xMax,Nquad) ; fQuad=xQuad+NaN ;
for k=1:Nquad
    fQuad(k)=J0+xQuad(k)*dJdp'*dp+xQuad(k)^2*dp'*Hess*dp/2 ;
end


x=linspace(0,xMax,N) ; f=x+NaN ;
parfor k=1:N
    f(k)=Func(x(k)) ;
end

f=f(:) ; x=x(:) ; fIn=fIn(:) ; xIn=xIn(:); 
f=[f;fIn]; x=[x;xIn]; 
[x,I]=sort(x) ; f=f(I); 


FindOrCreateFigure("TestQuadradicModel") ;
hold off
plot(xQuad,fQuad,'k--','LineWidth',2)
hold on
plot(x,f,'ro-')

legend('Quadradic approximation','function')

end