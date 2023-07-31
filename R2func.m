
function  R2=R2func(gamma,dx,dl,fun,x0,l0)

x=x0+gamma*dx;
l=l0+gamma*dl ;

R=fun(x) ;
R2=full(R'*R);



end

