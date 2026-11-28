function Blocks=AssembleHessianBlocks(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y,Type,K)

%%
%
%   Blocks=AssembleHessianBlocks(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y,Type)
%   Blocks=AssembleHessianBlocks(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y,Type,K)
%
% (16 Nov 2026) K: the forward-model Jacobian dF/dq at the state F, if available (eg from the forward solve in JGH.m, JGHouts.dFduv).
% It is then used directly, instead of being obtained from a further forward (uv) solve. F must then be the converged forward
% solution for which K was calculated, as returned by JGH.m.
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
%   Blocks.Kqq     : 2N x 2N, sparse part of K_qq
%   Blocks.KqqAG, Blocks.Kqqw : (16 Nov 2026) for dh/dt measurements in the -points- or -blocks- formulations, the dh/dt part of
%                    K_qq is (A G)' diag(w) (A G), with A G = [A Gu  A Gv] (number of measurements x 2N). It is kept in this factored
%                    form, and not formed: for block averages each row of A G contains all nodes of a block, so that the product
%                    couples all nodes within each block and (A G)' W (A G) has very many non-zeros (eg 1.8e8 for Antarctica nx150
%                    and 40 km blocks), while A G itself is sparse. Empty otherwise (then K_qq is all in Blocks.Kqq).
%   Blocks.Kpq     : Np x 2N
%   Blocks.Kpp     : Np x Np
%   Blocks.dA      : factorisation of K (MATLAB decomposition or PARDISO factorisation), reused for all solves in HessianVectorProduct.m
%   Blocks.dAT     : factorisation of K' (only if K is not symmetric)
%   Blocks.Reduced : (16 Nov 2026) if K is symmetric and the constraints are nodal: the free unknowns (Free) and the (Cholesky)
%                    factorisation of K(Free,Free) (Fac); then used instead of dA and dAT
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

narginchk(10,12)

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

%% forward-model Jacobian K=dF/dq at the current state, and the homogeneous constraints
% If not given, K is obtained from a forward (uv) solve starting from the state F (as in duv_hdABC.m).
if nargin<12 || isempty(K)
    [~,~,F,~,K]=uv([],[],CtrlVar,MUA,BCs,F,l);
end

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
KqqAG=[] ; Kqqw=[] ;
if is_dhdt_meas && dhdtMisfitMethod(CtrlVar)~="-nodal-integral-"
    % velocity part from Jqq (state independent, cached there), dh/dt part in factored form (see above)
    CtrlVarUV=CtrlVar ; CtrlVarUV.Inverse.Measurements=replace(CtrlVar.Inverse.Measurements,"-dhdt-","-") ;
    if contains(CtrlVarUV.Inverse.Measurements,"-uv-")
        Kqq=Jqq(CtrlVarUV,MUA,F,BCs,Meas) ;
    else
        Kqq=sparse(2*N,2*N) ;
    end
    [~,~,~,~,Mats]=dhdtObservationMisfit(CtrlVar,MUA,F,Meas,"Jacobians");
    KqqAG=[Mats.AGu Mats.AGv] ;
    Kqqw=Mats.w(:) ;
else
    Kqq=Jqq(CtrlVar,MUA,F,BCs,Meas) ;
end
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
Kpp=Jpp(CtrlVar,MUA,F,Meas,Type=="-GN-") ;   % (18 Nov 2026) Gauss-Newton: without the residual term of the dh/dt misfit (JhdotBB.m)
if Type=="-full-"
    Kpp=Kpp+Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
end

%% factorisation of K, reused for all solves
Blocks.Symmetric=norm(K-K',1)<=1e-12*norm(K,1) ;
frhs=zeros(2*N,1) ;
if isempty(L) ; grhs=[] ; else ; grhs=zeros(size(L,1),1) ; end
% (16 Nov 2026) The factorisation of K is returned by solveKApe (a MATLAB decomposition, or a PARDISO factorisation, see
% UaPardisoSolve.m) and reused in all solves of HessianVectorProduct.m, so that each Hessian-vector product only requires the solve
% phases. If K is not symmetric, K' is factorised as well, for the adjoint solves. Decompositions are supported by solveKApe for the
% constraint methods "Bempty" and "EliminateBCsSolveSystemDirectly"; for other methods dA=[] and each solve factorises the matrix.
%
% (16 Nov 2026) If K is symmetric (SSTREAM with Weertman sliding, as tested in SSTREAM2dNR2.m) and each constraint fixes one unknown
% (or there are no constraints), the same symmetric reduction as in the forward uv solve (solveKApeSymmetric.m, null-space method,
% KKT_null_space_decomposition_chol.m) is used: with the homogeneous constraints of the Hessian-vector products, x(c)=0 at the
% constrained unknowns c, and K(f,f) x(f) = r(f) at the free unknowns f. K(f,f) is symmetric (positive definite) and is factorised
% once, with a Cholesky factorisation (PARDISO "spd", which falls back to LDL' if K(f,f) is not positive definite, or MATLAB
% decomposition), and used for both the forward and the adjoint solves. Otherwise the general (nonsymmetric) route of solveKApe is used.
Blocks.Reduced=[] ; dA=[] ; dAT=[] ;
if Blocks.Symmetric
    [iRow,iCol,dVal]=find(L) ; m=size(L,1) ;
    if m==0 || (numel(dVal)==m && numel(unique(iRow))==m && numel(unique(iCol))==m)
        Free=true(2*N,1) ; Free(iCol)=false ;
        Kff=K(Free,Free) ; Kff=(Kff+Kff')/2 ;
        if isfield(CtrlVar,"LinearSolver") && isfield(CtrlVar.LinearSolver,"Type") && CtrlVar.LinearSolver.Type=="PARDISO"
            SPf=[] ;
            if isfield(MUA,"uvAssemblyPattern") && ~isempty(MUA.uvAssemblyPattern)
                SPf=struct("Pattern",MUA.uvAssemblyPattern,"Index",find(Free)) ;   % pattern of K(f,f), see UaPardisoSolve.m
            end
            Fac=UaPardisoSolve("factor",CtrlVar,Kff,[],"spd",SPf) ;
        else
            try
                Fac=decomposition(Kff,"chol") ;
            catch
                Fac=decomposition(Kff,"ldl") ;
            end
        end
        Blocks.Reduced=struct("Free",Free,"Fac",Fac) ;
    end
end

if isempty(Blocks.Reduced)
    [~,~,dA]=solveKApe(K,L,frhs,grhs,[],[],CtrlVar) ;
    if ~Blocks.Symmetric
        [~,~,dAT]=solveKApe(K',L,frhs,grhs,[],[],CtrlVar) ;
    end
end

Blocks.K=K ; Blocks.L=L ; Blocks.Fp=Fp ;
Blocks.Kqq=Kqq ; Blocks.KqqAG=KqqAG ; Blocks.Kqqw=Kqqw ; Blocks.Kpq=Kpq ; Blocks.Kpp=Kpp ;
% (23 Nov 2026) spatially correlated velocity errors: the correction -A K^{-1} A to the velocity part of K_qq (for u and v), applied
% in HessianVectorProduct.m, see CorrelatedPointDataMisfit.m (the same factors for u and v if the errors are equal)
Blocks.KqqCorrU=[] ; Blocks.KqqCorrV=[] ;
% (25 Nov 2026) optional prior on the area mean of log10(C): rank-one term w w'/sigma^2 in the C block (CLevelPrior.m)
[Blocks.CLevelW,Blocks.CLevelSigma]=CLevelPrior(CtrlVar,MUA) ;
Blocks.CLevelOffset=MUA.Nnodes*(contains(CtrlVar.Inverse.InvertFor,"logA")+contains(CtrlVar.Inverse.InvertFor,"-B-")) ;
if contains(CtrlVar.Inverse.Measurements,"-uv-") && is_uv_points_misfit(CtrlVar) && ~isempty(Meas.uvCorr)
    [~,~,Blocks.KqqCorrU]=CorrelatedPointDataMisfit(CtrlVar,MUA,Meas.uvO,Meas.usData,Meas.usErr,F.ub+F.ud,Meas.uvCorr) ;
    [~,~,Blocks.KqqCorrV]=CorrelatedPointDataMisfit(CtrlVar,MUA,Meas.uvO,Meas.vsData,Meas.vsErr,F.vb+F.vd,Meas.uvCorr) ;
end
Blocks.dA=dA ; Blocks.dAT=dAT ;
Blocks.Type=Type ; Blocks.Np=Np ; Blocks.N=N ;
Blocks.CtrlVar=CtrlVar ;

end
