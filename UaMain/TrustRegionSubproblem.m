

function alpha=TrustRegionSubproblem(H,E,g,alpha,Delta)

%%
% Solves a trust-region problem, i.e. finds a $\alpha$ such that 
% 
% $$\|p\| =  \Delta $$
%
% where
%
% $$ ( H + \alpha \, E ) p  = - g $$
%
% for given matrices $H$, $E$ and a vector $g$.
%
% $H$ and $E$ must be symmetrical matrices. 
%
% The above formulation is for $E$ being a diagonal matrix.
%
% If $E$ is not a diagonal matrix, the solution is to 
%
% $$ ( E^{-1} ( H + \alpha E ) ) \; p = - E^{-1}  g $$
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
% $H + \alpha E$ , and also needs to be guarded against the system becoming positive indefinite.
%
% The most relevant paper appears to be:  https://digital.library.unt.edu/ark:/67531/metadc283525/m2/1/high_res_d/metadc283525.pdf
% 
% But the method is also described in section 4.5 of Numerical Optimizations by Jorge Nocedal and Stephen J, Wright 
%
% https://www.ccom.ucsd.edu/~peg/papers/trust.pdf
%
%
%
% In the limiting case for $\alpha=0$ the system 
%
% $$ ( E^{-1} ( H + \alpha E ) ) \; p = - E^{-1}  g $$
% 
% becomes
%
% $$ H  p = - g $$
%
% which is the Newton system.
%
%
% And for $\alpha \to +\infty$  we have
%
%
% $$ \alpha \;  E   \; p = -  g $$
% 
% which give the direction of steepest descent. 
%
%
% The initial $\alpha$ value entered in the call may have to be a reasonable value for $\alpha$ for the method to converge.
%
% If the eigenvalues of $H$ are the set ${\lambda_i}$ where $\lambda_n \le \lambda_{n-1} \le \ldots \le \lambda_1$, 
% then $\alpha$ is in the interval $(-\lambda_n, +\infty)$.  Note that $H$ may not be positive definite and therefore
% we could have that $\lambda_1<0$.
%
% If $\alpha$ converge to a negative root in $(-\lambda_n, +\infty)$ the Newton step is inside the trust region and we can
% set $\alpha=0$ on return. This is, for example, the case when $\lambda_n>0$ and $H$ is positive definite and the Newton
% step inside the trust region. 
%
% There is an additional special case when the gradient $g$ is perpendicular to the eigenspace, $S_n$,  of the smallest $H$ eigenvalue, i.e. 
%
% $$ S_n= \{z : Hz =\lambda_n z \, , \, z \not = 0 \}$$
%
% This is sometimes referred to as the degenerate or the hard case.
%
% In the degenerate case, if $\lambda_n>0$ we can select $\alpha=0$ since $H$ is then positive definite and $\|p\|$, and if
% $\lamnbda_n$ is negative or zero, we can select $\alpha=\lambda_n \ge 0 $.
%
%
% For step size close to the full Newton step, the alpha values approach zero and Cholesky decomposition may not always work
% if the smallest eigenvalue of H are numerically small, yet positive. 
%
%
%
% Expansion on this using sub-space minimisation is the truncated Lanczos approach where $p$ is found as
%
% $$p_k=Q_k h_k$$ 
%
% see https://www.numerical.rl.ac.uk/media/people/nick-gould/GoulLuciRomaToin99_siopt.pdf
%
% Also see: https://people.maths.ox.ac.uk/nakatsukasa/preprints/TRSrev2.pdf
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



%% saveguarding 

lambdaMin=-1 ; % This should be the smallest eigenvalue of H, but this takes time to calculate...
alphaU=+inf;  % this should be smaller or equal to the smallest eigenvalue of H
alphaL=-inf ;

alphaLast=nan ;
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
%alpha=max(alpha,alphaL);

phiVector=nan(10,1);
alphaVector=nan(10,1);
phi=NaN; 

for iLoop=1:10

    HlE = H + alpha * E;

    [R, flag] = chol(HlE);  % change this later to use permutation matrix P once I switch to sparse Hessian


    if flag~=0
        fprintf("H+ alpha E not pos definite for alpha=%g \n",alpha)

        alphaL=alpha; % if not positive definite for this value of alpha, make sure not to use a lower value of alpha again.

        if isfinite(alphaU)
            alpha=max(0.001*alphaU,sqrt(alphaU*alphaL));
        else
            alpha=10*alphaL;
        end

        continue
    end

    p=R\(R'\(-g)) ;
    q=R'\p;

    pNorm=norm(p);
    qNorm=norm(q);


    phi=1/Delta-1/pNorm;

    phiVector(iLoop)=phi;
    alphaVector(iLoop)=alpha ; 

    if isEdiag
        dphidalpha=-(qNorm/pNorm)^2 * Delta/(pNorm-Delta) *phi;
    else
        dphidalpha=(phi-phiLast)/(alpha-alphaLast);
    end

    phiLast=phi; alphaLast=alpha;


    fprintf("===== it %i: phi=%-15g \t |p|=%-10g \t |q|=%-10g \t alpha=%g \n",iLoop,phi,pNorm,qNorm,alpha)



    if ~isnan(dphidalpha)
        alpha=alpha-phi/dphidalpha;
    else
        dalpha=1;
        alpha=alpha+dalpha;
    end


% the idea is to bracket the solution using the fact that phi is a monotonically decreasing function of alpha
%
% If phi >0 then I am to the left of the minimum, and I can set the lower bound for alpha, i.e. alphaL, to alphaLast; 
% If phi<0 then I am to the right of the minimum, and I can set the upper bound for alpha, i.e. alphaU, to alphaLast 
%

    if phi< 0 && alphaLast<alphaU
        alphaU=alphaLast;
    end


    if phi > 0 && alphaLast > alphaL % 
        alphaL=alphaLast;
    end

    if alpha > alphaU
        alpha=alphaU;
    end

    if alpha < alphaL
        alpha=alphaL;
    end

    % alpha=max(alpha,alphaL);
    % alpha=min(alpha,alphaU);
    % if alpha <= alphaS
    %     alpha=max(1e-6*alphaU,sqrt(alphaL*alphaU));
    % end

    if abs(phi)<phiTol
        fprintf("Tolerance for abs(phi) met. Exiting seach loop. \n")
        break
    end

end

alpha=alphaLast;
fprintf("Returning alpha=%g with phi=%g after %i iterations \n \n \n ",alpha,phi,iLoop)


%% plot phi

N=10;
alphaMax=1.5*alpha; 
alphaPlotVector=linspace(alphaL,alphaMax,N);
phiPlotVector=NaN(N,1);

for ILoop=1:N


    HlE = H + alphaPlotVector(ILoop) * E;

    [R, flag] = chol(HlE);  % change this later to use permutation matrix P


    if flag~=0
        fprintf("H+ alpha E not pos definite for alpha=%g.\n",alpha)
        continue
    end

    p=R\(R'\(-g)) ;

    pNorm=norm(p);

    phiPlotVector(ILoop)=1/Delta-1/pNorm;


end

FIG=FindOrCreateFigure("phi(alpha)") ; clf(FIG)
hold off
plot(alphaVector,phiVector,Color="r",Marker="o",MarkerFaceColor="b",LineStyle="none")
hold on

plot(alphaPlotVector,phiPlotVector,"or")


xlabel("$\alpha$",Interpreter="latex")
ylabel("$1/\Delta-1/\|p\|$",Interpreter="latex")
title("$ (H+\alpha E) p = -g $",Interpreter="latex")
yline(0,"--")

alphaPlotVector=[alphaPlotVector(:);alphaVector(:)];
phiPlotVector=[phiPlotVector(:);phiVector(:)];
I=~isnan(alphaPlotVector) | ~isnan(phiPlotVector);
alphaPlotVector=alphaPlotVector(I);
phiPlotVector=phiPlotVector(I);
[alphaPlotVector,I]=sort(alphaPlotVector);
phiPlotVector=phiPlotVector(I);

plot(alphaPlotVector,phiPlotVector,"-k")


plot(alpha,phi,Color="r",Marker="pentagram",MarkerFaceColor="b",MarkerEdgeColor="r",MarkerSize=22)

xline(alpha,"--","$\alpha^{\star}$",Interpreter="latex")

% XL=xlim;
% xline(alphaU,"--","$\alpha_U$",Interpreter="latex")
% xline(alphaL,"--","$\alpha_L$",Interpreter="latex")
% xlim(XL)

end