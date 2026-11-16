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
W=SolveK(Blocks.K,L,-Blocks.Fp*V,g,CtrlVar,Blocks.dA) ;

% 2)
A=Blocks.Kqq*W+Blocks.Kpq'*V ;
if isfield(Blocks,"KqqAG") && ~isempty(Blocks.KqqAG)
    A=A+Blocks.KqqAG'*(Blocks.Kqqw.*(Blocks.KqqAG*W)) ;   % dh/dt part of K_qq in factored form (16 Nov 2026)
end

% 3) adjoint solve
if Blocks.Symmetric
    Y=SolveK(Blocks.K,L,A,g,CtrlVar,Blocks.dA) ;
else
    Y=SolveK(Blocks.K',L,A,g,CtrlVar,Blocks.dAT) ;
end

% 4)
HV=-Blocks.Fp'*Y+Blocks.Kpq*W+Blocks.Kpp*V ;

end


function x=SolveK(K,L,f,g,CtrlVar,dA)

% solve K x = f subject to L x = g (several right-hand sides), reusing the decomposition dA if available

f=full(f) ;
if isempty(dA)
    x=solveKApe(K,L,f,g,[],[],CtrlVar) ;
else
    x=solveKApe(K,L,f,g,[],[],CtrlVar,dA) ;
end

end
