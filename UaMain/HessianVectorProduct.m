function HV=HessianVectorProduct(Blocks,V)

%%
%
%   HV=HessianVectorProduct(Blocks,V)
%
% Hessian-vector products H*V, for one or several vectors (the columns of V, Np x m), without forming the Hessian or the sensitivity
% matrix xi=dq/dp. Blocks is returned by AssembleHessianBlocks.m. With
%
% $$ H = \xi^T K_{qq} \xi + K_{pq} \xi + (K_{pq} \xi)^T + K_{pp} , \qquad \xi = -K^{-1} F_p $$
%
% the product is calculated as
%
%   1) linearised forward solve:   K W = -F_p V      (W = xi V), with homogeneous constraints L W = 0
%   2)                             A = K_qq W + K_pq' V
%   3) adjoint solve:              K' Y = A          (xi' A = -F_p' Y), with homogeneous constraints
%   4)                             H V = -F_p' Y + K_pq W + K_pp V
%
% ie two linear solves with K (and K'), plus sparse matrix-vector products. The factorisation of K (and of K' if K is not symmetric)
% is computed once in AssembleHessianBlocks.m and reused here, so that each product only requires the solve phases (16 Nov 2026).
%
% see also: AssembleHessianBlocks.m, CalcDirectAdjointHessian.m
%
%%

narginchk(2,2)

CtrlVar=Blocks.CtrlVar ;
m=size(V,2) ;
L=Blocks.L ;
if isempty(L) ; g=[] ; else ; g=zeros(size(L,1),m) ; end

% 1) W = xi V
W=SolveK(Blocks,Blocks.K,L,-Blocks.Fp*V,g,CtrlVar,Blocks.dA) ;

% 2)
A=Blocks.Kqq*W+Blocks.Kpq'*V ;
if isfield(Blocks,"KqqAG") && ~isempty(Blocks.KqqAG)
    A=A+Blocks.KqqAG'*(Blocks.Kqqw.*(Blocks.KqqAG*W)) ;   % dh/dt part of K_qq in factored form (16 Nov 2026)
end
if isfield(Blocks,"CLevelW") && ~isempty(Blocks.CLevelW)
    % (25 Nov 2026) rank-one Hessian of the prior on the area mean of log10(C), w w'/sigma^2, see CLevelPrior.m (added below to Hv)
    CLevelUsed=true;
else
    CLevelUsed=false;
end
if isfield(Blocks,"KqqCorrU") && ~isempty(Blocks.KqqCorrU)
    % (23 Nov 2026) spatially correlated velocity errors: -A K^{-1} A applied to the u and v parts, see CorrelatedPointDataMisfit.m
    Nn=size(Blocks.KqqCorrU.A,1) ;
    Cu=Blocks.KqqCorrU ; Cv=Blocks.KqqCorrV ;
    A(1:Nn,:)=A(1:Nn,:)-Cu.A*(Cu.L'\(Cu.L\(Cu.A*W(1:Nn,:)))) ;
    A(Nn+1:2*Nn,:)=A(Nn+1:2*Nn,:)-Cv.A*(Cv.L'\(Cv.L\(Cv.A*W(Nn+1:2*Nn,:)))) ;
end

% 3) adjoint solve
if Blocks.Symmetric
    Y=SolveK(Blocks,Blocks.K,L,A,g,CtrlVar,Blocks.dA) ;
else
    Y=SolveK(Blocks,Blocks.K',L,A,g,CtrlVar,Blocks.dAT) ;
end

% 4)
HV=-Blocks.Fp'*Y+Blocks.Kpq*W+Blocks.Kpp*V ;
if CLevelUsed
    iC=Blocks.CLevelOffset+(1:numel(Blocks.CLevelW)) ;
    HV(iC,:)=HV(iC,:)+Blocks.CLevelW*((Blocks.CLevelW'*V(iC,:))/Blocks.CLevelSigma^2) ;
end

end


function x=SolveK(Blocks,K,L,f,g,CtrlVar,dA)

% solve K x = f subject to the homogeneous constraints L x = 0 (several right-hand sides), with the symmetric reduced factorisation
% (Blocks.Reduced) if available, and otherwise reusing the decomposition dA if available

f=full(f) ;

if isfield(Blocks,"Reduced") && ~isempty(Blocks.Reduced)
    Free=Blocks.Reduced.Free ;
    x=zeros(size(f)) ;
    if isa(Blocks.Reduced.Fac,"decomposition")
        x(Free,:)=Blocks.Reduced.Fac\f(Free,:) ;
    else
        x(Free,:)=UaPardisoSolve("solve",Blocks.Reduced.Fac,f(Free,:)) ;
    end
    return
end

if isempty(dA)
    x=solveKApe(K,L,f,g,[],[],CtrlVar) ;
else
    x=solveKApe(K,L,f,g,[],[],CtrlVar,dA) ;
end

end
