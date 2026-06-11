function [x,y]=TestingMinres(A,B,f,g,x0,y0,CtrlVar)

% This is the simplest and most straigtforward iterative approach

n=size(A,1) ; m=size(B,1);
C=sparse(m,m);
AA=[A B' ;B -C] ;
bb=[f;g];

tluinc=tic;
P=cholinc(A,'0');  % the lower the droptolerance is the less number of iterations are needed
tluinc=toc(tluinc) ; disp([' cholinc in  ',num2str(tluinc),' sec ']) ;

R=[P spalloc(n,m,0) ; spalloc(m,n,0) speye(m,m)] ;



tol=1e-4 ; maxit=100;


[sol,flag,relres,iter,resvec]=minres(AA,bb,tol,maxit,R',R,[x0;y0]);

x=sol(1:n) ; y=sol(n+1:n+m);



semilogy(0:iter,resvec/norm(bb),'-o')
xlabel('Iteration Number')
ylabel('Relative Residual')

end



