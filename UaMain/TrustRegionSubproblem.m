

function lambda=TrustRegionSubproblem(H,E,g,lambda,Delta)

%%
% Solves a trust-region problem, i.e. finds a $\lambda$ such that 
% 
% $$\|p\| =  \Delta $$
%
% where
%
% $$ ( H + \lambda E ) p  = - g $$
%
% for given matrices $H$, $E$ and a vector $g$.
%
% $H$ and $E$ must be symmetrical matrices. 
%
% The above formulation is for $E$ being a diagonal matrix.
%
% If $E$ is not a diagonal matrix, the solution is to 
%
% $$ ( E^{-1} ( H + \lambda E ) ) \; p = - E^{-1}  g $$
% 
% with 
%
% $$\| p \| = \Delta $$
%
% This is equivalent to  solving the constraint quadratic minimization problem 
%
% $$ \min Q(p) = f + g' p + \frac{1}{2} p' H p $$
%
% subject to 
%
% $$ \| p \| \le \Delta $$
%
% The algorithm proves an almost exact solution of this problem, but does potentially require several (2 to 3) factorizations of
% $H + \lambda E$ , and also needs to be guarded against the system becoming positive indefinite.
%
% The most relevant paper appears to be:  https://digital.library.unt.edu/ark:/67531/metadc283525/m2/1/high_res_d/metadc283525.pdf
% 
% But the method is also described in section 4.5 of Numerical Optimizations by Jorge Nocedal and Stephen J, Wright 
%
%
% In the limiting case for $\lambda=0$ the system 
%
% $$ ( E^{-1} ( H + \lambda E ) ) \; p = - E^{-1}  g $$
% 
% becomes
%
% $$ H  p = - g $$
%
% which is the Newton system.
%
%
% And for $\lambda \to +\infty$  we have
%
%
% $$ \lambda \;  E   \; p = -  g $$
% 
% which give the direction of steepest descent. 
%
%
% The initial $\lambda$ value entered in the call may have to be a reasonable value for $\lambda$ for the method to converge.
%
%
%%

if isempty(E)

    nNodes=size(MUA.M,1);
    nH=size(H,1);

    if nH==nNodes
        E=MUA.M ;
    elseif nH==2*nNodes
        E=blkdiag(MUA.M,MUA.M) ;
    else
        error("wrong dimentions")
    end

end

% Can I use algorithm 4.3 to find lambda? Yes, but only if E is a diagonal matrix!
%
% https://digital.library.unt.edu/ark:/67531/metadc283525/m2/1/high_res_d/metadc283525.pdf
%

lambdaS=0 ; % lower bound on -lambda1 where lambda1 is the smallest eigenvalue of H. (this may be negative if H is not pos def.)
lambdaU=1e8 ;
lambdaL=0.01  ;

lambdaLast=nan ;
phiLast=nan;

n=numel(g);


isEdiag=isdiag(E);

if ~isEdiag

    H=E\H;
    g=E\g;
    H=0.5*(H+H');
    E=sparse(1:n,1:n,ones(n,1));
    isEdiag=true;
end


phiTol=1e-3; 
%lambda=max(lambda,lambdaL);

for iLoop=1:5

    HlE = H + lambda * E;

    [R, flag] = chol(HlE);  % change this later to use permutation matrix P once I switch to sparse Hessian


    if flag~=0
        fprintf("H+ l E not pos definite.\n")
    end

    p=R\(R'\(-g)) ;
    q=R'\p;

    pNorm=norm(p);
    qNorm=norm(q);


    phi=1/Delta-1/pNorm;

  

    if isEdiag
        dphidlambda=-(qNorm/pNorm)^2 * Delta/(pNorm-Delta) *phi;
    else
        dphidlambda=(phi-phiLast)/(lambda-lambdaLast);
    end

    phiLast=phi; lambdaLast=lambda;


    fprintf("===== it %i: phi=%-15g \t |p|=%-10g \t |q|=%-10g \t lambda=%g \n",iLoop,phi,pNorm,qNorm,lambda)



    if ~isnan(dphidlambda)
        lambda=lambda-phi/dphidlambda;
    else
        dlambda=1;
        lambda=lambda+dlambda;
    end


    % lambda=max(lambda,lambdaL);
    % lambda=min(lambda,lambdaU);
    % if lambda <= lambdaS
    %     lambda=max(1e-6*lambdaU,sqrt(lambdaL*lambdaU));
    % end

    if abs(phi)<phiTol
        break
    end

end

fprintf("Returning lambda=%g with phi=%g\n",lambda,phi)


%% plot phi

N=10;
lambdaMax=3*lambda; 
lambdaVector=linspace(lambda/100,lambdaMax,N);
phiVector=NaN(N,1);

for ILoop=1:N


    HlE = H + lambdaVector(ILoop) * E;

    [R, flag] = chol(HlE);  % change this later to use permutation matrix P


    if flag~=0
        fprintf("H+ l E not pos definite.\n")
        continue
    end

    p=R\(R'\(-g)) ;

    pNorm=norm(p);

    phiVector(ILoop)=1/Delta-1/pNorm;


end

FIG=FindOrCreateFigure("phi(lambda)") ; clf(FIG)

plot(lambdaVector,phiVector,"o-r")
xlabel("$\lambda$",Interpreter="latex")
ylabel("$1/\Delta-1/\|p\|$",Interpreter="latex")
title("$ (H+\lambda E) p = -g $",Interpreter="latex")
yline(0,"--")

hold on 
plot(lambda,phi,Color="r",Marker="o",MarkerFaceColor="b")

xline(lambda,"--","$\lambda^{\star}$",Interpreter="latex")



end