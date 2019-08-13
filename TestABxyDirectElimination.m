%%
n=5;
p=2;

A=rand(n,n) ; A=A+A' ; 
L=eye(p,n)  ;

f=ones(n,1) ; g=ones(p,1) ;


T=[A L' ; L zeros(p,p)];

sol=T\[f;g ] ; x=sol(1:n) ; y=sol(n+1:end) ; 


Q=(eye(n,n)-L'*L) ; 
xTest=(Q*A+ L'*L) \ (Q *f + L'*g)
(Q*A+ L'*L)

[x(:) xTest(:)]