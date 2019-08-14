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

%%
L=[1 0 0 0 ; 0 1 0 0 ] ; L=L./sqrt(sum(abs(L).^2,2)) ;isequal(L*L',eye(2,2))
L=[1 0 0 0 ; 0 0 1 0 ] ; L=L./sqrt(sum(abs(L).^2,2)) ;isequal(L*L',eye(2,2))
L=[1 0 0 0 ; 0 0 0 1 ] ; L=L./sqrt(sum(abs(L).^2,2)) ;isequal(L*L',eye(2,2))
L=[0 0 1 0 ; 0 1 0 0 ] ; L=L./sqrt(sum(abs(L).^2,2)) ;isequal(L*L',eye(2,2))
L=[0 0 0 1 ; 0 1 0 0 ] ; L=L./sqrt(sum(abs(L).^2,2)) ;isequal(L*L',eye(2,2))

%%

%%

L=[1 1 0 -1 ; 0 0 1 0 ]  ; L=L./sqrt(sum(abs(L).^2,2)) ; norm(L*L'-eye(2,2)) < 100*eps

%%
L=[ 1 10 0 -2 0 0 0  ; ...
    0 0 2  0 0  0 0  ; ...
    0 0 0  0 5  0 0  ; 
    0 0 0  0 0  1 -1]; 

L=L./sqrt(sum(abs(L).^2,2)); 
L*L'

p=size(L,1) ; 
norm(L*L'-eye(p,p)) < 100*eps

%%
syms x
L=[ 1 x 0 -2 0 0 0  ; ...
    0 0 2  0 0  0 0  ; ...
    0 0 0  0 5  0 0  ; 
    0 0 0  0 0  1 -1]; 
