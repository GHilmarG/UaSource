
function  R2=R2func(gamma,dx,fun,x0)


narginchk(4,4)

x=x0+gamma*dx;
%l=lambda0+gamma*dlambda ;

R=fun(x) ;
R2=full(R'*R);



end

