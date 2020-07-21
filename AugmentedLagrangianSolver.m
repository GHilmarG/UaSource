function [x,y] = AugmentedLagrangianSolver(A,B,f,g,y0,CtrlVar)


isUpperLeftBlockMatrixSymmetrical=issymmetric(A);


if CtrlVar.InfoLevelLinSolve>=10
    if isUpperLeftBlockMatrixSymmetrical
        fprintf(' Solving a symmetrical indefinite block system using the Augmented Lagrangian Solver (ALS) \n' )
    else
        fprintf(' Solving asymetrical indefinite block system using the Augmented Lagrangian Solver (ALS) \n' )
    end
end

IterationMin=CtrlVar.ALSIterationMin;
IterationMax=CtrlVar.ALSIterationMax;

% solves a system on the form [A B' ; B 0]=[f;g]
% using Augmented-Lagrangian method, where the inner problem is solved directly
%
% LU/LDL factorisation done outside of loop, hence cost of each additional iteration fairly small
% The inner iteration is so cheap in comparision to the LU/LDL factorisation
% that it seems justified to always do at least two iterations
%
% Method:
% There seems to be confusion in the literature about what the Uzawa and the ALS methods are.
% For
% [A B'] [x]= [f]
% [B 0 ] [y]  [g]
% the Uzawa method is somtimes defined as:
%                       A x_{i+1}= f - B' y_i
%                         y_{i+1}=g y_i + iW B x_{i+1}
%  where iW is `small' but not too small.
% 
%  I find that the following modified version of the Uzawa methods converges very quickly:
%
%  Define T= [A B']
%            [B iW]
% where iW is `small'.
%
%  Repeat: T [x_{i+1}] = [    f    ]
%            [y_{i+1]]   [g+ iW y_i]
% 
% This idea of replacing the lower-left 0 block matrix with something small, and redefining the
% rhs accordingly, is usually referred to as the Augmented-Lagrangian method. 
% Note that the iteration, once converged, gives 
% answer that is a correct solution to the unmodifed system, and that this final answer does not
% depend on what iW is chosen to be. So iW does not need to be `small'. However, usually 
% the convergence is fastest if iW is small compared to A.  If setting iW=0 works,
% then only one solve is needed (fastest possible rate of convergence).
%
%  One of the nice aspects of the iteration is that T has only to be factorized once.
% (LDL if symmetric, LU otherwise) and each additional iteration is cheap


%[m,n]=size(B);
n=size(A,1) ; m=size(B,1);
x0=zeros(n,1);



% The following seems a good way of selecting iW
%k=round(log10(mean(abs(diag(A))))) ; w=10^(k+ALSpower);  % w=1;
k=round(log10(norm(diag(A))))  ;   
w=10^(k+CtrlVar.ALSpower)  ;
iW=speye(m)/w;



T=[A B' ; B iW];

tStart=tic;

luvector=CtrlVar.Solve.LUvector;


if isUpperLeftBlockMatrixSymmetrical &&  CtrlVar.TestForRealValues
    [L,D,p,S]=ldl(T,'vector');   % LDL factorisation using MA57, MA57 is a multifrontal sparse direct solver using AMD ordering
    sol=zeros(m+n,1);
elseif luvector
    
    if isdistributed(T)
        T=full(T);                 % can only be done for full martices 
        [L,U,p] = lu(T,'vector');  % apparently this can be used in a parallel mode, but used in non-parallel this is very slow!
        % sol=U\(L(p,:)\fg) ;
    else
        [L,U,p,q,R] = lu(T,'vector');  % this can not be used in a parallel mode
        % sol(q)=U\(L\(R(:,p)\fg)) ;
    end
    
    sol=zeros(m+n,1);
    
else
    
    if isdistributed(T)
        [L,U,P] = lu(T); % lu factorisation not using UMFPACK, much slower!
        % T xy = fg and  T=P' L U - > xy=U\(L\(P*fg)) 
        % anyhow, this is not implemented for sparse distributed arrays anyhow
    else
        [L,U,P,Q,R] = lu(T); % lu factorisation using UMFPACK
    end
    
end


% requires 2k^2 n operations where n is the matrix size and k the semi-bandwidth
tElapsed=toc(tStart);
if CtrlVar.InfoLevelLinSolve>1 ; fprintf(' LU factorisation in %g sec \n',tElapsed) ; end

%     lu factorisation for the non-symmetrical case
%     Q*Q'=1,  P*P'=1,  but  R*R' ~=1
%     A = R*P'*L*U*Q'
%     P*(R\A)*Q = L*U
%     sol=Q*(U\(L\(P*(R\fg))));
%
% LDL factorisation for the symmetrical case
% P'*S*A*S*P = L*D*L'
% T=S\(P*L*D*L'*P')/S)
% S*A*S = P*L*D*L'*P'
% S=S' , P*P'=1
% sol=S*P*(L'\(D\(L\(P'*S*fg))));
%
% LU vector
%  Example: 
%
%   A x = y ; 
%  [L,U,p,q,R] = lu(A,'vector');
%  x(q)=U\(L\(R(:,p)\y)) 
%
%    T sol = fg
%    [L,U,p,q,R] = lu(T,'vector');
%    sol(q)=U\(L\(R(:,p)\fg)) 


% I measure the residual as res=norm([A B' ; B 0]-[f;g])/norm([f ; g])
% It is possible that the norm(B*x-g)/norm(y) will become comparable to
% machine precision and then the iteration should stop, hence the need for res2

Iteration=0;

if CtrlVar.InfoLevelLinSolve>=10 
    InfoVector=zeros(IterationMax+1,5) ;
end

resRelative=1e10 ; xDiff=1e10; yDiff=1e10; resAbsolute=1e10;

if isempty(y0) ; y0=zeros(m,1) ; end

while (resRelative > CtrlVar.LinSolveTol &&  resAbsolute > 1e-10 && Iteration <= IterationMax) || Iteration < IterationMin
    % I continue as long as either relative or absolute residual is too large
    % Even if the relative residual is not smaller, it makes no sense to continue if the absolute residual is approaching numerical resolution
    
    Iteration=Iteration+1;
    fg=[f ; g + iW*y0];
    % sol=Q*(U\(L\(P*(R\fg))));
    
    if isUpperLeftBlockMatrixSymmetrical &&  CtrlVar.TestForRealValues
        fg=S*fg ; sol(p)=L'\(D\(L\(fg(p)))); sol=S*sol;  % if using the vector format
        
    elseif luvector
        if isdistributed(T)
            sol=U\(L(p,:)\fg) ; % not sure if correct, could not test because lu not yet implemented for distributed sparse matrices!
        else
            sol(q)=U\(L\(R(:,p)\fg)) ;
        end
    else
        if isdistributed(T)
            sol=U\(L\(P*fg)) ;  % not sure if correct, could not test because lu not yet implemented for distributed sparse matrices!
        else
            sol=Q*(U\(L\(P*(R\fg))));   % P*(R\A)*Q = L*U for sparse non-empty A.
        end
    end
    
    if any(isnan(sol))
        save TestSave
        error('AugmentedLagrangianSolver:NaN','NaN in sol. All variables writen to TestSave.mat')
    end
    
    x=sol(1:n) ; y=sol(n+1:end);
    
    %[A B'] [x]=[f]
    %[B 0 ] [y]=[g]
    resAbsolute=norm([A B' ; B sparse(m,m)]*sol-[f ; g])/norm([f;g]);
    resRelative=resAbsolute/norm([f;g]);
    
    
    yDiff=norm(y-y0)/norm(y);
    xDiff=norm(x-x0)/norm(x);
    %res=norm([A B' ; B sparse(m,m)]*sol-[f ; g]);
    
    if CtrlVar.InfoLevelLinSolve>=10 
        res1=norm(A*x+B'*y-f)/norm(f);
        ng=norm(g);
        if ng>0
            res2=norm(B*x-g)/norm(g);
        else
            res2=norm(B*x-g);
        end
        
        
        InfoVector(Iteration,1)=xDiff;
        InfoVector(Iteration,2)=yDiff;
        InfoVector(Iteration,3)=resRelative;
        InfoVector(Iteration,4)=res1;
        InfoVector(Iteration,5)=res2;
    end
    x0=x;
    y0=y;
    
    if xDiff<eps && yDiff< eps 
        % sometimes the solution no longer changes
        % this could be a real issue, but it might also numerical resolution issue
        % I here force the normalized residuals tolerance to be no smaller than 1e-8, which seems a
        % very save assumption. 
        CtrlVar.LinSolveTol=max([CtrlVar.LinSolveTol ; 1e-8]) ; 
        break
    end
    
end

if CtrlVar.InfoLevelLinSolve>=10
    fprintf(' Number of Augmented Lagrangian Iterations when solving ([A B'' ; B 0] [x;y] =[f;g]) was  %-i \n',Iteration);
    fprintf(' Relative solution residual %-g \t , relative change in x and y: %-g and %-g \n',resRelative,xDiff,yDiff);
end

if Iteration > IterationMax
    res1=norm(A*x+B'*y-f);
    res2=norm(B*x-g);
    if CtrlVar.InfoLevelLinSolve<1
        fprintf(' Relative and absolute total solution residuals:  %-g , %-g  \n',resRelative,resAbsolute);
        fprintf('  Change in x and y: %-g and %-g \n',xDiff,yDiff);
    end
    fprintf(' Absolute solution residuals for first equation %-g \t and second equation %-g \n',res1,res2);
    res1=res1/norm(f);
    res2=res2/norm(g);
    fprintf(' Relative solution residuals for first equation %-g \t and second equation %-g \n',res1,res2);
    warning('ALS:MaxIterationReached','Augmented Lagrangian Solver exits because maximum number of iterations %g reached \n',IterationMax)
end

if resRelative > CtrlVar.LinSolveTol &&  resAbsolute > 1e-10
    res1=norm(A*x+B'*y-f)/norm(f);
    
    ng=norm(g);
    if ng>0
        res2=norm(B*x-g)/norm(g);
    else
        res2=norm(B*x-g);
    end
    fprintf(' relative residuals=%-g \t absolute residuals=%-g \t first equation %-g \t second equation %-g \n',resRelative,resAbsolute,res1,res2);
    warning('ALS:MaxIterationReached','Augmented Lagrangian Solver did not fully converge to prescribed tolerance of %-g \n',CtrlVar.LinSolveTol)
    
    % save TestSave
    
end

if CtrlVar.InfoLevelLinSolve>=10
    fprintf(' --------------------------- AugmentedLagrangianLinSolver: ----------------------------------------------------------------------------------------------\n ')
    fprintf(' solves a system on the form [A B'' ; B 0] [x;y] =[f;g]  iterativily \n')
    fprintf(' starting with [A B'' ; B iW]=[f;g]  , where iW=1/10^(k+CtrlVar.ALSpower)  using k=%g and CtrlVar.ALSpower=%g \n ',...
        full(k),CtrlVar.ALSpower)
    ng=norm(g);
    if ng>0
        fprintf(' Iteration           |x-x0|/|x|                  |y-y0|/|y|           |[A B'' ; B 0]*sol-[f;g]|/|[f;g]|    |A*x+B''*y-f|/|f|        |B*x-g|/|g| \n')
    else
        fprintf(' Iteration           |x-x0|/|x|                  |y-y0|/|y|           |[A B'' ; B 0]*sol-[f;g]|/|[f;g]|    |A*x+B''*y-f|/|f|        |B*x-g|\n')
    end
    for I=1:Iteration
        fprintf('%10i \t %20.10g \t \t %20.10g \t \t %30.10g \t %20.10g \t %20.10g \n',...
            I,InfoVector(I,1),InfoVector(I,2),InfoVector(I,3),InfoVector(I,4),InfoVector(I,5))
    end
    fprintf(' --------------------------------------------------------------------------------------------------------------------------------------------------------\n ')
    iRange=1:Iteration;  
    figure ; 
    semilogy(iRange,InfoVector(iRange,1),'+-r'); 
    hold on ; 
    semilogy(iRange,InfoVector(iRange,2),'+-g'); 
    semilogy(iRange,InfoVector(iRange,3),'o-k'); 
    semilogy(iRange,InfoVector(iRange,4),'o-m'); 
    semilogy(iRange,InfoVector(iRange,5),'o-c'); 
    legend('|x-x0|/|x|','|y-y0|/|y|','|[A B'' ; B 0]*sol-[f;g]|/|[f;g]|','|A*x+B''*y-f|/|f|','|B*x-g|/|g|') ; 
    title(' Augmented-Lagrangian linear solver ([A B'' ; B 0] [x;y] =[f;g])')
    xlabel(' Iteration '); ylabel(' ')
    input(' to continue press return ')   

end


end



