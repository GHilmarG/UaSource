

function KFpp=Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

%% Builds the Hessian term
%
% $$ \mathcal{F}^{pp} = \langle \Psi , \delta^2_{pp} \mathcal{F} \rangle $$
%
% i.e. the second-order derivatives of the forward model with respect to the inversion fields
% $p=(\log_{10} A, B, \log_{10} C)$, contracted with the adjoint variables.
%
%% Block structure
%
% Which blocks are non-zero follows directly from where each field enters the forward model:
%
%  * the viscous terms contain $\eta(A)$ and the thickness $h$, and $h$ depends on $B$ through the geometrical
%    closure. Hence $\mathcal{F}^{AB} \neq 0$.
%  * the basal drag contains $\beta^2(C)$ and the grounding mask $\mathcal{G}$, and $\mathcal{G}$ depends on $B$.
%    (For effective-pressure dependent sliding laws there is a further route through $N$.) Hence
%    $\mathcal{F}^{BC} \neq 0$.
%  * no term in the forward model contains both $A$ and $C$: the viscosity involves no $C$, and no sliding law
%    involves $A$. Hence $\mathcal{F}^{AC} = 0$ exactly.
%
% The full matrix is therefore block tridiagonal in the ordering $(A,B,C)$
%
% $$ \mathcal{F}^{pp} = \left[\begin{array}{ccc}
%  \mathcal{F}^{AA} & \mathcal{F}^{AB} & 0                \\
%  \mathcal{F}^{BA} & \mathcal{F}^{BB} & \mathcal{F}^{BC} \\
%  0                & \mathcal{F}^{CB} & \mathcal{F}^{CC}
% \end{array}\right] $$
%
% with $\mathcal{F}^{BA}=(\mathcal{F}^{AB})^T$ and $\mathcal{F}^{CB}=(\mathcal{F}^{BC})^T$.
%
% B sits in the middle because it enters through the geometry, which multiplies both the viscous and the basal-drag
% terms. When B is not inverted for, the matrix reduces to the block-diagonal form used previously.
%
% The blocks corresponding to inactive fields are dropped at the end, so any combination of the three fields is
% handled by the same code.
%
%  see also: FAA.m, FBB.m, FCC.m, FAB.m, FBC.m, Jpp.m, CalcDirectAdjointHessian.m
%
%%

narginchk(7,7)
nargoutchk(1,1)

[isA,isB,isC]=isABC(CtrlVar);

Active=[isA isB isC] ;

if ~any(Active)
    error("Fpp:NoActiveFields","None of the inversion fields A, B or C is active.")
end

%% The cross blocks involving B
%
% The A-B and B-C blocks are calculated by FAB.m and FBC.m, respectively, see below. Both require CtrlVar.Hh0=0 and no melange model physics,
% and FBC.m is only implemented for the Weertman sliding law (as FCuv.m and Fqq.m).

%% assemble the blocks

nNodes=MUA.Nnodes ;

K=repmat({sparse(nNodes,nNodes)},3,3) ;   % all blocks default to zero

% diagonal blocks

if isA
    K{1,1}=FAA(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

if isB
    K{2,2}=FBB(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

if isC
    K{3,3}=FCC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

% off-diagonal blocks
%
% Note: K{1,3} and K{3,1}, i.e. the A-C blocks, are identically zero and are therefore left at their default value.

if isA && isB
    K{1,2}=FAB(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;  K{2,1}=K{1,2}.' ;
end

if isB && isC
    K{2,3}=FBC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;  K{3,2}=K{2,3}.' ;
end

%% keep only the active fields

KFpp=cell2mat(K(Active,Active)) ;

end
