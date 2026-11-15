function Blocks=AssembleHessianBlocks(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y,Type)

%%
%
%   Blocks=AssembleHessianBlocks(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y,Type)
%
% Assembles the sparse matrices needed to calculate Hessian-vector products (see HessianVectorProduct.m), without forming the
% (dense) sensitivity matrix xi=dq/dp or the (dense) Hessian itself.
%
% The Hessian of the cost function with respect to the parameters p (logA, B, logC), with q=(u,v) the velocities and F(q,p)=0 the
% forward model, is (see CalcDirectAdjointHessian.m)
%
% $$ H = \xi^T K_{qq} \xi + K_{pq} \xi + (K_{pq} \xi)^T + K_{pp} , \qquad \xi = -K^{-1} F_p $$
%
% where K=dF/dq, F_p=dF/dp, and
%
%   Type="-GN-"   (Gauss-Newton, default):  K_qq = J^{qq} ,              K_pq = J^{pq}_GN ,          K_pp = J^{pp}
%   Type="-full-" :                         K_qq = J^{qq} + F^{qq}[Psi] , K_pq = J^{pq} + F^{pq}[Psi] , K_pp = J^{pp} + F^{pp}[Psi]
%
% Here J^{pq}_GN is the Gauss-Newton part of the mixed derivative of the dh/dt misfit term (non-zero only for B inversions with dh/dt
% measurements, see JhdotBq.m), and Psi the adjoint variables. These are the same terms as calculated by CalcDirectAdjointHessian.m
% with HessianTerms="-xi Jqq xi-Jpp-" (GN) and "-all-" (full), so that H*v from HessianVectorProduct.m agrees with those Hessians.
%
% The Gauss-Newton Hessian is positive semi-definite. The full Hessian is in general indefinite away from the minimum (and can be
% slightly indefinite even at the minimum).
%
% On return, Blocks contains:
%
%   Blocks.K       : 2N x 2N  forward-model Jacobian dF/dq
%   Blocks.L       : constraint matrix (velocity boundary conditions), the HVP uses the homogeneous constraints L x = 0
%   Blocks.Fp      : 2N x Np  dF/dp (with respect to log10(A), B, log10(C))
%   Blocks.Kqq     : 2N x 2N
%   Blocks.Kpq     : Np x 2N
%   Blocks.Kpp     : Np x Np
%   Blocks.dA      : decomposition of K (if available), reused for all linear solves in HessianVectorProduct.m
%   Blocks.Symmetric : true if K is symmetric (then the adjoint solves reuse the decomposition of K)
%   Blocks.Type, Blocks.Np, Blocks.N
%
% The blocks depend on the current state (F, and Psi for Type="-full-"), and must be reassembled when p changes.
%
% Limitations as for CalcDirectAdjointHessian.m: SSTREAM, Weertman sliding law for the C terms, no box transform.
%
% see also: HessianVectorProduct.m, CalcDirectAdjointHessian.m, duv_hdABC.m, Hess_qp.m
%
%%

narginchk(10,11)

if nargin<11 || isempty(Type)
    Type="-GN-";
end
Type=string(Type);
if ~ismember(Type,["-GN-","-full-"])
    error("AssembleHessianBlocks:Type","Type=%s, but must be either -GN- or -full-.",Type)
end

if CtrlVar.Inverse.BoxTransform
    error("AssembleHessianBlocks:NotImplemented","The box transform is not implemented for the Hessian blocks (as for CalcDirectAdjointHessian.m).")
end

[isA,isB,isC] = isABC(CtrlVar);
[~,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar) ;
N=MUA.Nnodes ;
Np=(isA+isB+isC)*N ;
ln10=log(10) ;

%% forward-model Jacobian K=dF/dq at the current state (as in duv_hdABC.m), and the homogeneous constraints
[~,~,F,~,K]=uv([],[],CtrlVar,MUA,BCs,F,l);

BCs0=BCs ;
if numel(BCs0.ubFixedValue)>0 ; BCs0.ubFixedValue=BCs0.ubFixedValue*0 ; end
if numel(BCs0.vbFixedValue)>0 ; BCs0.vbFixedValue=BCs0.vbFixedValue*0 ; end
L=AssembleLuvSSTREAM(CtrlVar,MUA,BCs0) ;

%% F_p = dF/dp, blocks ordered (A,B,C), with respect to log10(A) and log10(C) (as in duv_hdABC.m)
FpA=[] ; FpB=[] ; FpC=[] ;
if isA ; FpA=dFuvdA(CtrlVar,MUA,F)*spdiags(F.AGlen(:)*ln10,0,N,N) ; end
if isB ; FpB=dFuvdBGeneral(CtrlVar,MUA,F) ; end
if isC ; FpC=dFuvdC(CtrlVar,MUA,F)*spdiags(F.C(:)*ln10,0,N,N) ; end
Fp=[FpA FpB FpC] ;

%% K_qq
Kqq=Jqq(CtrlVar,MUA,F,BCs,Meas) ;
if Type=="-full-"
    Kqq=Kqq+Fqq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

%% K_pq (Np x 2N), rows ordered (A,B,C), columns (u,v), as in Hess_qp.m
if Type=="-full-"
    [KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
    [KFBu,KFBv]=FBuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
    [KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
    Kpq=[KFAu KFAv ; KFBu KFBv ; KFCu KFCv] ;
else
    Kpq=sparse(Np,2*N) ;
end
if isB && is_dhdt_meas
    iB=isA*N+(1:N) ;
    Kpq(iB,:)=Kpq(iB,:)+JhdotBq(CtrlVar,MUA,F,Meas,Type=="-GN-") ;
end

%% K_pp
Kpp=Jpp(CtrlVar,MUA,F,Meas) ;
if Type=="-full-"
    Kpp=Kpp+Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

%% factorisation of K, reused for all solves
Blocks.Symmetric=norm(K-K',1)<=1e-12*norm(K,1) ;
frhs=zeros(2*N,1) ;
if isempty(L) ; grhs=[] ; else ; grhs=zeros(size(L,1),1) ; end
dA=[] ;
try
    [~,~,dA]=solveKApe(K,L,frhs,grhs,[],[],CtrlVar) ;   % returns the decomposition used, if supported (not with PARDISO)
catch
    dA=[] ;
end

Blocks.K=K ; Blocks.L=L ; Blocks.Fp=Fp ;
Blocks.Kqq=Kqq ; Blocks.Kpq=Kpq ; Blocks.Kpp=Kpp ;
Blocks.dA=dA ;
Blocks.Type=Type ; Blocks.Np=Np ; Blocks.N=N ;
Blocks.CtrlVar=CtrlVar ;

end
