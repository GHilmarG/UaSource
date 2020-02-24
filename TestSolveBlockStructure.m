%%
%%
fprintf('\n\n\n')
n=3 ;
p=2;

A=sprandsym(n,n,0.01);
A=sparse([1 0.01 0 ; 0.01 1 0 ; 0 0 1]);


B=zeros(p,n) ;

for I=1:p
    B(I,I)=1;
end

C=zeros(p,p);
AA=[A B' ; B C ]  ;

f=zeros(n,1)+1;
g=zeros(p,1)+1;

% Null-space method
tic
Z=null(B);
xhat=B\g ;
z= (Z'*A*Z)\(Z'*(f-A*xhat));
x=Z*z+xhat;
y=B' \ (f-A*x) ;
toc

fprintf('%g \t %g \n ',norm(A*x+B'*y-f),norm(B*x-g))

% Range-space method
tic
y= B*(A\B')\(B*(A\f)-g) ;
x=A\(f-B'*y) ;
toc



%  R=chol(A)        % A = R'*R    ->  A x = b ; 
%                                  % x = A\b 
%                                  % x = R\(R'\b)
%                                  % A\b = (R\R'\b)
% 
% R=chol(A) ;                          
%  [R,flag,p]= chol (A, 'vector') ->  A(p,p) = R'*R  
%  [R,flag,P]= chol (A )          -> P'*A*P = R'*R  -> A\b = (
%                            
%
%
% [R,flag,P]= chol (A); 
% sol=A\f ; 
% sol2=P*(R\(R'\(P'*f))) ; 
% norm(sol-sol2)


%% this works for chol and sparse
 tic
 
 % y= B*(R\(R'\B'))\(B*(R\(R'\f))-g) ;
 % x= R\(R'\(f-B'*y)) ;
 [R,flag,P]= chol (A); 
 y= B*P*(R\(R'\(P'*B')))\(B*P*(R\(R'\(P'*f)))-g) ;
 x= P*(R\(R'\(P'*f-P'*B'*y))) ;
 
 toc                              

fprintf('%g \t %g \n ',norm(A*x+B'*y-f),norm(B*x-g))
%%


tic
CtrlVar.ALSIterationMin=1;
CtrlVar.ALSIterationMax=10;
CtrlVar.ALSpower=5; CtrlVar.Solve.LUvector=true; 
CtrlVar.InfoLevelLinSolve=1;
CtrlVar.TestForRealValues=true;
CtrlVar.LinSolveTol=1e-6;
[x,y] = AugmentedLagrangianSolver(A,B,f,g,g*0,CtrlVar);
toc
fprintf('%g \t %g \n ',norm(A*x+B'*y-f),norm(B*x-g))

