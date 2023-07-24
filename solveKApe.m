function  [x,y]=solveKApe(A,B,f,g,x0,y0,CtrlVar)

narginchk(7,7)

% Solves:
%
%  [A   B'] [x]= [f]
%  [B  -C ] [y]  [g]
%
% where A is n times n, C is m times m , B is m times n.
% A does not have to be symmetrical

[nA,mA]=size(A) ; [nB,mB]=size(B) ; [nf,mf]=size(f) ; [ng,mb]=size(g) ;
[nx0,mx0]=size(x0) ; [ny0,my0]=size(y0);

if nA~=mA
    fprintf(' A must be square ')
end

if mB~=0 && (nA~=mB || mA ~= mB)
    fprintf('size of A (%-i,%-i) and B (%-i,%-i) matrices not consistent \n',nA,mA,nB,mB)
    error('error in solveKApe')
end

if nf~=nA
    fprintf('f must have same number of elements as there are rows in A \n')
    error('error in solveKApe')
end

if ng~=nB
    fprintf('g has %g but must have same number of elements as there are rows in B ie %g \n ',ng,nB)
    error('error in solveKApe')
end

if ~isempty(B)
    if nx0~=nA
        fprintf('x0 must have same number of elements as there are rows in A\n')
        error('error in solveKApe')
    end
    
    if ny0~=nB
        fprintf('y0 must have same number of elements as there are rows in B\n')
        error('solveKApe:InputsIncompatable','error in solveKApe')
    end
end

n=size(A,1) ; m=size(B,1);

%
% if isempty(B) || numel(B)==0
%     CtrlVar.AsymmSolver='Bempty';
% elseif all(full(sum(B~=0,2))==1)
%     %isequal(B*B',sparse(1:m,1:m,1))  % if only one node is constrained in each constraint, then pre-eliminate and solve directly
%     CtrlVar.AsymmSolver='EliminateBCsSolveSystemDirectly';
% end

if isempty(CtrlVar) || ~isstruct(CtrlVar)

    CtrlVar.AsymmSolver='auto';
    CtrlVar.InfoLevelLinSolve=0 ;
    CtrlVar.TestForRealValues=0;

else
    if ~isfield(CtrlVar,"AsymmSolver")
        CtrlVar.AsymmSolver='auto';
    end

    if ~isfield(CtrlVar,"InfoLevelLinSolve")
        CtrlVar.InfoLevelLinSolve=0 ;
    end

    CtrlVar.TestForRealValues=0;
end

if isequal(lower(CtrlVar.AsymmSolver),'auto')

    if isempty(B) || numel(B)==0
        CtrlVar.AsymmSolver='Bempty';
    elseif isdiag(B*B')
        CtrlVar.AsymmSolver='EliminateBCsSolveSystemDirectly';
    else
        CtrlVar.AsymmSolver='AugmentedLagrangian';
    end

end



tSolve=tic;

switch CtrlVar.AsymmSolver
    case 'Bempty'
        
        x=A\f; y=[];
        
    case 'Backslash'
        
        m=size(B,1); C=sparse(m,m); AA=[A B' ;B -C] ; bb=[f;g];
        sol=AA\bb; x=sol(1:n) ; y=sol(n+1:end);
        if CtrlVar.InfoLevelLinSolve>=1
            fprintf(' Constraint matrix NOT empty. Solving system directly using the backslash operator \n ')
        end
        
    case 'EliminateBCsSolveSystemDirectly'
    
        
         [x,y]=ABfgPreEliminate(CtrlVar,A,B,f,g);
    
    case 'AugmentedLagrangian'
        
        
        [x,y] = AugmentedLagrangianSolver(A,B,f,g,y0,CtrlVar);
        
    case 'EliminateBCsSolveSystemIterativly'
        
        if CtrlVar.InfoLevelLinSolve>2; fprintf(' Eliminating constraints and solving system iterativly \n') ; end
        
        [I,iConstrainedDOF]=ind2sub(size(B),find(B==1)); iConstrainedDOF=iConstrainedDOF(:);
        iFreeDOF=setdiff(1:n,iConstrainedDOF); iFreeDOF=iFreeDOF(:);
        
        
        AA=A; ff=f; xx0=x0;
        AA(iConstrainedDOF,:)=[]; AA(:,iConstrainedDOF)=[]; ff(iConstrainedDOF)=[];
        xx0(iConstrainedDOF)=[];
        
        tstart=tic;


        tluinc=tic;
        %setup.type = 'crout'; setup.milu = 'off'; setup.droptol = 0.1;
        %setup.type = 'ilutp'; setup.milu = 'off'; setup.droptol = 0.15;
        setup.type = 'nofill'; setup.milu = 'off';
        
        [L1,U1] = ilu(AA,setup);
        tluinc=toc(tluinc);
        
        
        tol=1e-6 ; maxit=20;
        


        t1=tic ;
        %[sol,flag,relres,iter,resvec]=bicgstabl(AA,ff,tol,maxit,L1,U1,xx0);
        restart=10;
        [sol,flag,relres,iter,resvec]=gmres(AA,ff,restart,tol,maxit,L1,U1,xx0);
        t2=toc(t1);
        
        %sol=AA\ff;
        x=zeros(n,1) ; x(iConstrainedDOF)=g ; x(iFreeDOF)=sol; y=B*(f-A*x);
        
        
        tend=toc(tstart);
        
        if CtrlVar.InfoLevelLinSolve>=1
            disp([' ilu in  ',num2str(tluinc),' sec '])
            disp([' gmres  ',num2str(t2),' sec '])
            disp([' total solution time  ',num2str(tend),' sec '])
        end
        
        if CtrlVar.InfoLevelLinSolve>=10
            
            
            figure
            fprintf(' flag=%-i, iter=%-g, relres=%-g \n ',flag,iter,relres)
            nnn=numel(resvec);
            semilogy((0:nnn-1)/2,resvec,'-o')
            xlabel('Iteration Number')
            ylabel('Relative Residual')
        end
        
    otherwise
        
        error(' which case ? ')
        
end

tSolve=toc(tSolve);
if CtrlVar.InfoLevelLinSolve>=10
    fprintf('solveKApe: # unknowns=%-i \t # variables=%-i \t # Lagrange mult=%-i \t time=%-g \t method=%s \n ',...
        n+m,n,m,tSolve,CtrlVar.AsymmSolver)
    if CtrlVar.InfoLevelCPU
        fprintf(CtrlVar.fidlog,' in %-g sec. \n',tSolve) ;
    end
end

%% Testing
if isempty(B)
    if norm(f)>eps
        res=norm(A*x-f)/norm(f);
    else
        res=norm(A*x-f);
    end
    if res>1e-5
        fprintf('solveKApe: Solution residual appears too large! %g \n',res)
    end
else
    if norm(f)>1e-10
        res1=norm(A*x+B'*y-f)/norm(B'*y-f);
    elseif norm(f)==0  && (norm(x) > 0 || norm(y)>0) 
        res1=norm(A*x+B'*y)/norm(B'*y);
    else
        res1=norm(A*x+B'*y-f);
    end
    if norm(g)>1e-10
        res2=norm(B*x-g)/norm(g);
    else
        res2=norm(B*x-g);
    end
    if res1> 1e-5 || res2 > 1e-5
        fprintf('solveKApe: Solution residuals appear to be too large! %g %g \n',res1,res2)
    end
end

if CtrlVar.TestForRealValues
    
    if ~isreal(x) && CtrlVar.IgnoreComplexPart
        x=real(x);
    end
    
    if ~isreal(y) && CtrlVar.IgnoreComplexPart
        y=real(y) ;
    end
    
end

% 		%[ Schur complement reduction
% 		%  -S:=B/A*B'+C  (negative Schur compliment)
% 		% (5.1):  (B/A*B'+C)*y=B/A*f-g;
% 		% (5.2):   A x=f-B'*y
% 		% When A and -S are symmetric positive definite, solve using chol or CG.
%
% 		C=sparse(m,m);
%
% 		%C=eps*speye(m);  % just a slight regularistation is sufficient to make S pos. definite
%
% 		if Schur==1
% 			tSchur=tic;
% 			T=B/A; % T=B*inv(A) = B/A  = (A'\B')'  so T=(A'\B')'=(A\B')' if A=A'
% 			S=T*B'+C ;
% 			y=S\(T*f-g);
% 			x=A\(f-B'*y);
%
%
% 			tSchur=toc(tSchur) ;
%
% 			AA=[A B' ;B -C] ; bb=[f;g];
% 			test=AA*[x;y]-bb; disp([' direct Schur : ',num2str(max(abs(test))),' in ',num2str(tSchur)])
%
% 		end
%
%
%
% 		%
% 		if AugmentedLagrangian==1
%
%
%
% 		end
%
%
%
% 		if Bicgstab==1   % seems to work fine with using ilu with setup.type='nofill';
% 			tstart=tic ;
% 			% no modification of system,
% 			AA=[A B' ;B -C] ; bb=[f;g];
% 			if mod(iteration,10)==1
%
% 				AAA=[A B' ;B  sparse(1:m,1:m,1)] ;
%
%
% 				tluinc=tic;
% 				%setup.type = 'crout'; setup.milu = 'row'; setup.droptol = 0.1;
% 				setup.type = 'nofill'; setup.milu = 'off';
%
% 				[L1,U1] = ilu(AAA,setup);
% 				%[L1,U1] = luinc(AAA,0);
% 				tluinc=toc(tluinc);
% 				disp([' ilu in  ',num2str(tluinc),' sec '])
% 			end
%
%
% 			tol=1e-7 ; maxit=50;
%
% 			if nargin<5 ; x0=zeros(n,1) ; y0=zeros(m,1) ; end
%
% 			%M1=sparse(1:n+m,1:n+m,[spdiags(A,0);ones(m,1)],n+m,n+m);
% 			%sol=bicgstab(AA,bb,tol,maxit,M1,[],[x0;y0]);
%
%
% 			t1=tic ;
%
% 			[sol,flag,relres,iter,resvec]=bicgstabl(AA,bb,tol,maxit,L1,U1,[x0;y0]);
%
%
% 			t2=toc(t1);
% 			disp([' bicgstab  ',num2str(t2),' sec '])
%
%
%
% 			x=sol(1:n) ; y=sol(n+1:end);
% 			tend=toc(tstart);
% 			disp([' total solution time  ',num2str(tend),' sec '])
%
% 			figure
% 			fprintf(' flag=%-i, iter=%-g, relres=%-g \n ',flag,iter,relres)
% 			nnn=numel(resvec);
% 			semilogy((0:nnn-1)/2,resvec/norm(bb),'-o')
% 			xlabel('Iteration Number')
% 			ylabel('Relative Residual')
%
% 		end
%
% 	end
% end
%
