%%
n=4;
p=2;

A=rand(n,n) ; A=A+A' ; 
B=eye(p,n)  ;

f=ones(n,1) ; g=ones(p,1) ;


T=[A B' ; B zeros(p,p)];

sol=T\[f;g ] ; x=sol(1:n) ; y=sol(n+1:end) ; 


Q=(eye(n,n)-B'*B) ; 
xTest=(Q*A+ B'*B) \ (Q *f + B'*g)

format long g
[x(:) xTest(:)]