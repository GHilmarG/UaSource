%%
%%
fprintf('\n\n\n')
n=350 ;
p=200;

A=sprandsym(n,n,0.01);

L=zeros(p,n) ;

for I=1:p
    L(I,I)=1;
end

B=zeros(p,p);
AA=[A L' ; L B ]  ;

b=zeros(n,1)+1;
c=zeros(p,1)+1;

% Null-space method
tic
Z=null(L);
xhat=L\c ;
z= (Z'*A*Z)\(Z'*(b-A*xhat));
x=Z*z+xhat;
lambda=L' \ (b-A*x) ;
toc

fprintf('%g \t %g \n ',norm(A*x+L'*lambda-b),norm(L*x-c))

% Range-space method
tic
lambda= L*(A\L')\(L*(A\b)-c) ;
x=A\(b-L'*lambda) ;
toc


% tic
% iA=inverse(factorize(A)) ;
% lambda= L*(iA*L')\(L*(iA*b)-c) ;
% x=iA*(b-L'*lambda) ;
% toc

fprintf('%g \t %g \n ',norm(A*x+L'*lambda-b),norm(L*x-c))

tic
CtrlVar.ALSIterationMin=1;
CtrlVar.ALSIterationMax=10;
CtrlVar.ALSpower=5; CtrlVar.Solve.LUvector=true; 
CtrlVar.InfoLevelLinSolve=1;
CtrlVar.TestForRealValues=true;
CtrlVar.LinSolveTol=1e-6;
[x,lambda] = AugmentedLagrangianSolver(A,L,b,c,c*0,CtrlVar);
toc
fprintf('%g \t %g \n ',norm(A*x+L'*lambda-b),norm(L*x-c))