function  [x,y]=solveKApeSymmetric(A,B,f,g,x0,y0,CtrlVar)

% Solves:
%
%  [A   B'] [x]= [f]
%  [B   0 ] [y]  [g]
%
% where A is nxn, C is mxm , B is mxn
% The system must be symmetric

%nA=size(A,1) ; nB=size(B,1);

% if isempty(y0) && ~isempty(B);  % if y0 is given as empty, set the initial estimate to zero
%     y0=zeros(size(B,1),1);
% end

[nA,mA]=size(A) ; [nB,mB]=size(B) ; [nf,mf]=size(f) ; [ng,mb]=size(g) ;  
%[nx0,mx0]=size(x0) ; 

if isempty(y0)
  y0=zeros(nB,1); % This is a special case, allowing for the initial estimate for y to be empty
                  % in which case the initial estimate is set to zero.
end

[ny0,my0]=size(y0);

if nA~=mA
    fprintf(' A must be square ')
end

if mB~=0 && (nA~=mB || mA ~= mB)
    fprintf('size of A (%-i,%-i) and B (%-i,%-i) matrices not consistent \n',nA,mA,nB,mB)
    save TestSave ; error('error in solveKApeSymmetric')
end

if nf~=nA
    fprintf('f must have same number of elements as there are rows in A \n')
    save TestSave ; error('error in solveKApeSymmetric')
end

if ng~=nB
    fprintf('g has %g but must have same number of elements as there are rows in B ie %g \n ',ng,nB)
    save TestSave ; error('error in solveKApeSymmetric')
end

% x0 is never used
% if nx0~=nA
%     fprintf('x0 must have same number of elements as there are rows in A\n')
%     error('error in solveKApeSymmetric')
% end



if ny0~=nB
    fprintf('y0 must have same number of elements as there are rows in B\n')
    save TestSave ; error('error in solveKApeSymmetric')
end

if ~isfield(CtrlVar,"Symmsolver")
   CtrlVar.SymmSolver='auto' ;
end


if isequal(lower(CtrlVar.SymmSolver),'auto')
    
    if isempty(B) || numel(B)==0
        CtrlVar.SymmSolver='Bempty';
    elseif isdiag(B*B')
        CtrlVar.SymmSolver='EliminateBCsSolveSystemDirectly';
    else
        CtrlVar.SymmSolver='AugmentedLagrangian';
    end
    
end

tSolve=tic; 

switch CtrlVar.SymmSolver
    case 'Bempty'
        x=A\f;
        y=[];
        
    case 'AugmentedLagrangian'

   
        CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=1;
        [x,y] = AugmentedLagrangianSolver(A,B,f,g,y0,CtrlVar);

        
    case 'Backslash'
        
        C=sparse(nB,nB);
        AA=[A B' ;B -C] ; bb=[f;g];
        sol=AA\bb;
        x=sol(1:nA) ; y=sol(nA+1:nA+nB);
        
    case 'EliminateBCsSolveSystemDirectly'
     
        [x,y]=ABfgPreEliminate(CtrlVar,A,B,f,g);

    otherwise
        
        error('case not reckognised ')
        
end


tSolve=toc(tSolve);

if isfield(CtrlVar,"InfoLevelLinSolve")
    if CtrlVar.InfoLevelLinSolve>=10
        fprintf('solveKApeSymmetric: # unknowns=%-i \t # variables=%-i \t # Lagrange mult=%-i \t time=%-g \t method=%s \n ',...
            nA+nB,nA,nB,tSolve,CtrlVar.SymmSolver)
    end
end


return

end


%
%
%
% if solutionphase==1 ;
%     x0=x0one ; y0=y0one ;
%     if InfoLevel==1 ; disp(' using x0one and y0one' ) ; end
% elseif solutionphase==2 ;
%     x0=x0two ; y0=y0two ;
%     if InfoLevel==1 ; disp(' using x0two and y0two' ) ; end
% end
%
% if isempty(x0)  ; x0=zeros(n,1) ; end
% if isempty(y0)  ; y0=zeros(m,1) ; end
%
% if ~all(x0==0) && ~all(y0==0)
%
%     res=norm(A*x0+B'*y0-f)/norm(x0) + norm(B*x0-g)/norm(x0);
%     tol=1e-10;
%     if res < tol ; disp([' solution norm ',num2str(res),' < ', num2str(tol),' solveKApeSymmetric returns ']) ;
%         x=x0 ; y=y0;
%         return
%     end
% end
%
%
%
% C=sparse(m,m);
%
%
%
% %  Remark about iteration: The system can be written as
% % x=(A-B' C B)\(f-B'y-B'C g)
% % y=C\(B x -g)
% % which can be used in a iterative fashion (I think) as
%
% % x_i=(A-B' C B)\(f-B'y_i-B'C g)
% % y_{i+1}=C\(B x_i -g)
% %
% % Now if I have a system on the form
% %  [A   B'] [x]= [f]
% %  [B   0 ] [y]  [g]
% % I can modify it and write
% %  [A   B'] [x]= [f]
% %  [B   -C] [y]  [g-C y]
% % which leads to an iterative approach
% %
% % x_i=(A-B' C B)\(f-B'y_i-B'C g)
% % y_{i+1}=y_i+C\(B x_i -g)
% % This is the iterative augmented Lagrangian method where C=I/w and w is `large'
%
%
% % currently best results with Uzawa=2 and Minres=3
%
%
%
% %     if all(x0==0) && all(y0==0)
% %         Minres=0;
% %         Schur=0;
% %         Uzawa=2;  % Uzawa=2 seems to give best results, no singularity (Augmented Lagrangian)
% %         SymmetricLQ=0;
% %         if InfoLevel>1 ; disp(' using Uzawa 2 ') ; end
% %     else
% %         %Minres=1;  % minres where the A matrix is modifed, but
% %         %Minres=2; % minres used with augmented Lagrangian stabilisation with outer Uzawi iteration (pos def)
% %
% %         Minres=0;  %=3 minres used on unmodified system
% %         Schur=0;
% %         Uzawa=2;  % Uzawa=2 seems to give best results, no singularity (Augmented Lagrangian)
% %         SymmetricLQ=0;
% %         if InfoLevel>1 ; disp(' using Uzawa 2') ; end
% %
% %     end
%
%
%
%
%
% %[ Schur complement reduction
% %  -S:=B/A*B'+C  (negative Schur compliment)
% % (5.1):  (B/A*B'+C)*y=B/A*f-g;
% % (5.2):   A x=f-B'*y
% % When A and -S are symmetric positive definite, solve using chol or CG.
%
%
%
% %C=eps*speye(m);  % just a slight regularistation is sufficient to make S pos. definite
%
% if strcmp(SolMethod,'Schur1')
%     tSchur=tic;
%     T=B/A; % T=B*inv(A) = B/A  = (A'\B')'  so T=(A'\B')'=(A\B')' if A=A'
%     S=T*B'+C ;
%     y=S\(T*f-g);
%     x=A\(f-B'*y);
%
%
%     tSchur=toc(tSchur) ;
%
%     AA=[A B' ;B -C] ; bb=[f;g];
%     test=AA*[x;y]-bb; disp([' direct Schur : ',num2str(max(abs(test))),' in ',num2str(tSchur)])
%
% end
%
%
%
% %[-------- Uzawa method
% if strcmp(SolMethod,'Uzawa1')
%
%     y1=zeros(m,1);
%
%     tUzawa=tic;
%     w=1/eigs(B/A*B',1)/2;
%     %w=1e6;  % 0 < w < 1/lambda_max  of B/A*B'
%     tol=1e-5; res=tol*10; IImax=2000; II=0;
%
%     while res > tol && II <= IImax
%
%         for I=1:10
%             II=II+1;
%             % Augmented Lagrangian
%             x1=(A-w*(B'*B))\(f-B'*y1-w*B'*g);
%             y1=y1+w*(B*x1-g);
%
%             %x1=A\(f-B'*y1);
%             %y1=y1+w*(B*x1-g);
%
%         end
%         res=max(abs(A*x1+B'*y1-f));
%         if InfoLevel==1 ; disp([' Uzawa  residual : ',num2str(res)]) ; end
%     end
%     if InfoLevel==1 ; disp([' Number of Uzawa iterations  ',num2str(II)]) ; end
%
%     tUzawa=toc(tUzawa);
%     if InfoLevel==1 ; disp([' Uzawa  in ',num2str(tUzawa)]) ; end
%     x=x1 ; y=y1;
% end
%
% if strcmp(SolMethod,'Minres1')
%     %AA=[A B' ;B -C] ; bb=[f;g];
%     k=round(log10(max(diag(A)))) ;
%     gamma=10^(k+3) ;
%     %gamma=1e8;
%     tol=1e-6 ; maxit=1000;
%     % only exact if C=0
%     Anew=A+gamma*(B'*B);
%     AA=[Anew B' ;B -C] ; bb=[f+gamma*B'*g;g]; % Augmented Lagrangian
%
%
%     M1=sparse(1:n+m,1:n+m,[spdiags(Anew,0);ones(m,1)],n+m,n+m);
%     %M1=[1e-6*A sparse(n,m); sparse(m,n) sparse(1:m,1:m,1e-6*ones(m,1),m,m)];
%
%
%     tminres=tic ; sol=minres(AA,bb,tol,maxit,M1,[],[x0;y0]); tminres=toc(tminres);
%     x2=sol(1:n) ; y2=sol(n+1:end);
%     if InfoLevel==1 ; disp([' minres solves in  ',num2str(tminres),' sec ']) ; end
%
%     x=x2 ; y=y2;
% end
%
%
% if strcmp(SolMethod,'SymmLQ')
%     %AA=[A B' ;B -C] ; bb=[f;g];
%     gamma=1e6;
%     AA=[A+gamma*(B'*B) B' ;B -C] ; bb=[f+gamma*B'*g;g]; % Augmented Lagrangian
%
%     tol=1e-8 ; maxit=1000;
%     tLQ=tic ; sol=symmlq(AA,bb,tol,maxit,[],[],[]); tLQ=toc(tLQ);
%     x2=sol(1:n) ; y2=sol(n+1:end);
%     if InfoLevel==1 ; disp([' Symmetric LW solves in  ',num2str(tLQ),' sec ']) ; end
%
%     x=x2 ; y=y2;
% end
%
%
%
%
% if strcmp(SolMethod,'Uzawa2')
%
%     UzawaIterationMin=3;
%     tol=1e-20; % tolerance etc for Uzawa iteration
%     [x,y] = UzawaSymmSolver(A,B,f,g,y0,tol,InfoLevel,UzawaIterationMin);
%
%
%
% end
%
%
% if strcmp(SolMethod,'Uzawa3')
%     % Uzawa with conjugated dirctions as outer interation (C=0)
%     % I used direct solver to solve a system with A
%     % should try a version were I use conjugated gradients to solve that one
%
%
%     tUzawa=tic;
%     tol=1e-10; IImax=30; II=0;
%
%
%     y=y0;
%     x=A\(f-B'*y) ;
%
%     d=B*x-g;
%     q=-d;
%
%     res=norm(q)/sqrt(length(q));
%     res2=norm(A*x+B'*y-f)/sqrt(length(f));
%     disp([' Uzawa3  residual : ',num2str(res),' ',num2str(res2)])
%
%     while res > tol  && II <= IImax
%
%         II=II+1;
%         p=B'*d;
%         h=A\p;
%         alpha=q'*q/(p'*h);
%         y=y+alpha*d;
%         x=x-alpha*h;
%         qold=q;
%         q=g-B*x;
%         beta=q'*q/(qold'*qold);
%
%         d=-q+beta*d;
%         res=norm(q)/sqrt(length(q));
%         res2=norm(A*x+B'*y-f)/sqrt(length(f));
%         disp([' Uzawa3  residual : ',num2str(res),' ',num2str(res2)])
%     end
%     tUzawa=toc(tUzawa);
%     disp([' Number of Uzawa3 iterations  ',num2str(II)])
%     disp([' Uzawa  solves in ',num2str(tUzawa),' sec '])
% end
%
% %%
%
% if strcmp(SolMethod,'Uzawa4')
%     % Uzawa with conjugated dirctions as outer iteration (C=0)
%     % and iterative inner solution as well
%     % seems that the outer it, requires the inner to be solve very accuratly
%     tUzawa=tic;
%
%     tolinner=1e-5; % this is reduced in the loop
%     tolouter=1e-5;
%
%     IImax=20; maxit=1000 ; II=0;
%
%
%     y=y0; x=x0;
%
%     if mod(iteration,50)==1 || isempty(R)
%         tcinc=tic;
%         R = cholinc(A,'0');
%         %P=cholinc(A,0.001);
%         tcinc=toc(tcinc); disp([' cholinc in  ',num2str(tcinc),' sec '])
%     end
%
%     tpcg=tic;
%     x = pcg(A,f-B'*y,tolinner,maxit,R',R,x); %x=A\(f-B'*y) ;
%     tpcg=toc(tpcg); disp([' pcg in  ',num2str(tpcg),' sec '])
%
%     d=B*x-g;
%     q=-d;
%
%     res=norm(q)/sqrt(length(q));
%     res2=norm(A*x+B'*y-f)/sqrt(length(f));
%     disp([' Uzawa4  residual : ',num2str(res),' ',num2str(res2)])
%
%     h=x;
%     while res > tolouter  && II <= IImax || II <1
%
%         II=II+1;
%         p=B'*d;
%         %h=A\p;
%         tolinner=tolinner/2;
%         h = pcg(A,p,tolinner,maxit,R',R,h);
%         alpha=q'*q/(p'*h);
%         y=y+alpha*d;
%         x=x-alpha*h;
%         qold=q;
%         q=g-B*x;
%         beta=q'*q/(qold'*qold);
%
%         d=-q+beta*d;
%         res=norm(q)/norm(x);
%         res2=norm(A*x+B'*y-f)/norm(x);
%         disp([' Uzawa4  residual : ',num2str(res),' ',num2str(res2)])
%     end
%     tUzawa=toc(tUzawa);
%     disp([' Number of Uzawa4 iterations  ',num2str(II)])
%     disp([' Uzawa  solves in ',num2str(tUzawa),' sec '])
% end
%
%
%
%
%
% if strcmp(SolMethod,'Minres2')
%     if InfoLevel>1 ; fprintf('Minres 2') ; end
%     % here I use the minres in combination with the iterative augmented Lagrangian
%     % because the matrix is  positive definite I can use the incomplete Cholesky
%     % factorisation for preconditioning
%
%     k=round(log10(max(diag(A)))) ; w=10^(k+6);  iW=speye(m)/w;
%     T=[A B' ; B iW];
%
%     if InfoLevel>1 ; fprintf('cholinc') ; end
%     R = cholinc(T,1e-3);
%     tol=1e-5; tolU=1e-5 ; IImax=4; II=0; maxit=1000; res=tolU*10;
%
%     while res > tolU  && II <= IImax
%
%
%         II=II+1;
%         fg=[f ; g + iW*y0];
%         %M1=sparse(1:n+m,1:n+m,[spdiags(T,0)],n+m,n+m);
%         %M1=sparse(1:n+m,1:n+m,[spdiags(A,0);ones(m,1)],n+m,n+m);
%         %tminres=tic ; sol=minres(T,fg,tol,maxit,M1,[],[x0;y0]); tminres=toc(tminres);
%
%         tminres=tic ; sol=minres(T,fg,tol,maxit,R',R,[x0;y0]); tminres=toc(tminres);
%         tol=tol/10;  % good to have lower tolerance in next iteration so that it does not just return the previous solution
%         x=sol(1:n) ; y=sol(n+1:end);
%         disp([' minres solves in  ',num2str(tminres),' sec '])
%
%         res=max(abs(y-y0));
%         if InfoLevel>1 ; disp([' Uzawa2  residual : ',num2str(res)]) ; end
%         y0=y; x0=x;
%
%     end
%     if InfoLevel>1 ; disp([' Number of outer Uzawa2 iterations  ',num2str(II)]) ; end
%
% end
%
% if strcmp(SolMethod,'Minres3')
%     % This is the simplest and most straigtforward iterative approach
%     iminres=iminres+1;
%     ttminres=tic;
%     % just use minres on the unmodified system!
%
%     if InfoLevel>1 ; disp(' Minres on unmodified system ') ; end
%
%     % no modification of system,
%     C=sparse(m,m);
%     AA=[A B' ;B -C] ;
%     bb=[f;g];
%
%
%     % approx to Schur complement: B (A\B')+C
%
%
%     if isempty(P) || mod(iminres,50)==0
%         if InfoLevel>1 ; fprintf('cholinc \n') ; end
%         tluinc=tic;
%         P=cholinc(A,'0');  % the lower the droptolerance is the less number of iterations are needed
%         tluinc=toc(tluinc); if InfoLevel>1 ; disp([' cholinc in  ',num2str(tluinc),' sec ']) ; end
%     end
%
%     %T=speye(m,m);  % I tried to approximate the Schur complemet, but that was slower!
%     R=[P spalloc(n,m,0) ; spalloc(m,n,0) speye(m,m)] ;
%
%     % R=speye(n+m,n+m); % no preconditioning
%
%
%
%
%     tol=1e-8 ; maxit=n;
%     tminres=tic ;
%     if InfoLevel>1 ; fprintf('minres \n') ; end
%     [sol,flag,relres,iter,resvec]=minres(AA,bb,tol,maxit,R',R,[x0;y0]);
%     tminres=toc(tminres);
%     if InfoLevel>1 ; disp([' minres solves in  ',num2str(tminres),' sec ']) ; end
%     if flag~=0 ; disp([ ' minres flag : ',num2str(flag)]) ; endfunction  [x,y]=solveKApeSymmetric(A,B,f,g,x0,y0,CtrlVar)
%
