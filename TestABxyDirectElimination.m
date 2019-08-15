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
%%

%load solveKApeTestSave     ; % just velocity contraints, works
Klear
% load solveKApeTestSaveWithThicknessConstraints
%load solveKApeTestSaveWithManyThicknessConstraints
load solveKApeTestSavePIGtransient

% new pre-elimination followed by a direct solve
p=size(B,1) ; Test=B*B'-speye(p,p); sum(sum(abs(Test))) ; 
n=size(A,1) ;

Q=speye(n,n)-B'*B;
Atilde=Q*A+B'*B;
btilde=(Q*f+B'*g) ;
tic
for I=1:5
    xTest=Atilde\btilde; 
    yTest=B*(f-A*xTest);
end
preEliminationTime = toc ;

tic
for I=1:5
    dAtilde=factorize(Atilde);
    xTest3=dAtilde\btilde;
    yTest3=B*(f-A*xTest3);
end
preEliminationAutoFactorizeTime = toc ;

tic
for I=1:5
    [LL,UU,PP,QQ,RR] = lu(Atilde); 
    xTest2=QQ*(UU\(LL\(PP*(RR\btilde))));
    yTest2=B*(f-A*xTest2);
end
preEliminationManualFactorizeTime = toc ;

% matlab backslash
p=size(B,1); C=sparse(p,p); AA=[A B' ;B -C] ; bb=[f;g];
tic
for I=1:5
    sol=AA\bb;
    xBack=sol(1:n) ; yBack=sol(n+1:end);
end
BackslashTime = toc ;
%
tic
for I=1:5
    [xAL,yAL] = AugmentedLagrangianSolver(A,B,f,g,y0,CtrlVar);
end
AugmentedLagrangianTime = toc ;
%

fprintf('\n preElimination \t \t \t \t %f \n preEliminationAutoFactorize \t %f \n preEliminationManualFactorize \t %f \n Backslash \t \t \t \t \t \t %f \n AugmentedLagrangian \t \t \t %f \n',...
    preEliminationTime,preEliminationAutoFactorizeTime,preEliminationManualFactorizeTime,BackslashTime,AugmentedLagrangianTime)

[norm(xTest-xBack) norm(xTest2-xBack) norm(xTest3-xBack) norm(xAL-xBack)]



