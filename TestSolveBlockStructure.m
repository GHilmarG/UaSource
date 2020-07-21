%%
%%
fprintf('\n\n\n')
nA=5 ;
p=2;


A=sparse([1 0.01 0 0 0.1  ; 0.01 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0.1 ; 0.1 0 0 0.1 1]);

% singular A  (some of the methods still work even if A does not have info provided B provides that
% info)
% A=sparse([0 0.0 0 0 0  ; 0 0 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0.1 ; 0 0 0 0 1]);


B=zeros(p,nA) ;

for I=1:p
    B(I,I)=1;
end

C=zeros(p,p);
AA=[A B' ; B C ]  ;

f=zeros(nA,1)+1;
g=zeros(p,1)+1;

%%
fprintf('\n\n\n')

[nA,mA]=size(A); 

% Null-space method  (does not work for sparse
if ~issparse(B)
    tic
    Z=null(B);
    xhat=B\g ;
    z= (Z'*A*Z)\(Z'*(f-A*xhat));
    x=Z*z+xhat;
    y=B' \ (f-A*x) ;
    fprintf(' Null-space:  \t \t \t \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))
end

% Range-space method
tic
y= (B*(A\B'))\(B*(A\f)-g) ;
x=A\(f-B'*y) ;


fprintf('Range-space:  \t \t \t \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))


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


% this works for chol and sparse
tic

% y= B*(R\(R'\B'))\(B*(R\(R'\f))-g) ;
% x= R\(R'\(f-B'*y)) ;
[R,flag,P]= chol (A);
y= B*P*(R\(R'\(P'*B')))\(B*P*(R\(R'\(P'*f)))-g) ;
x= P*(R\(R'\(P'*f-P'*B'*y))) ;



fprintf('Range-space (chol):  \t \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))

%

% Range-space method if A singular
tic
A=A+B'*B ; 
f=f+B'*g;
y= (B*(A\B'))\(B*(A\f)-g) ;
x=A\(f-B'*y) ;


fprintf('Range-space (singular): \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))

%%

%%
%
% Pre-elimination, not symmetric
for I=1:1
    tic
    if isdiag(B*B')
        
        factor=1./(full(mean(abs(diag(A)))));
        A=factor*A ; f=factor*f ;  % this leaves x unaffected but y is scaled
        % factor=1;
        Q=speye(nA,nA)-B'*B;
        
        
        Atilde=Q*A+B'*B;
        btilde=(Q*f+B'*g) ;

        
        %dAtilde=factorize(Atilde);
        
        x=Atilde\btilde;
        y=B*(f-A*x);
        
        A=A/factor ; f=f/factor ;  % this leaves x unaffected but y is scaled
        y=y/factor;
        
        
        fprintf('Pre-eliminate:  \t \t \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))
    end
end

%

% ALS


tic
CtrlVar.ALSIterationMin=1;
CtrlVar.ALSIterationMax=10;
CtrlVar.ALSpower=5; CtrlVar.Solve.LUvector=true;
CtrlVar.InfoLevelLinSolve=1;
CtrlVar.TestForRealValues=true;
CtrlVar.LinSolveTol=1e-6;
[x,y] = AugmentedLagrangianSolver(A,B,f,g,g*0,CtrlVar);
fprintf('ALS:  \t \t \t \t \t \t %f sec \t %g \t %g \n ',toc,norm(A*x+B'*y-f)/norm(f),norm(B*x-g))

