






function [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC,Meas)

narginchk(13,14)
nargoutchk(1,1)

%% Builds the mixed Hessian term
%
% $$ H^{pq}+H^{qp} = \big( \mathcal{F}^{pq} \big) \xi + \Big[ \big( \mathcal{F}^{pq} \big) \xi \Big]^T $$
%
% With velocity measurements only, the $J^{pq}$ contribution is absent because each term in the cost function is an explicit function of either p or
% q, but not of both.
%
% If dh/dt measurements are included, $J_{\dot{h}}$ depends on the velocities AND on the thickness, and hence on B (through h=h(B)), so $J^{pq} \neq 0$
% for B. This term is added to the B row of $\mathcal{F}^{pq}$ below, see JhdotBq.m, and Meas must then be provided as the (optional) last input.
% For A and C inversions (ie where p=A or p=C) the corresponding contribution is identically equal to zero, because they do not affect h.
% 
%
% The blocks are ordered (A,B,C), matching Fpp.m and the assembly of xi in CalcDirectAdjointHessian.m
%
% Inactive fields drop out automatically: FAuv.m, FBuv.m and FCuv.m each return empty matrices when their field is
% not being inverted for, and the corresponding sensitivity matrices are empty as well, so the concatenations below
% simply omit those rows and columns.
%
%  see also: FAuv.m, FBuv.m, FCuv.m, Fpp.m, CalcDirectAdjointHessian.m
%
%%

[isA,isB,isC] = isABC(CtrlVar);

[~,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar) ;
if isB && is_dhdt_meas && (nargin<14 || isempty(Meas))
    error("Hess_qp:MeasNeeded","Hess_qp needs the input Meas when inverting for B with dh/dt measurements.")
end


[KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

[KFBu,KFBv]=FBuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

[KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;


KFpq=[KFAu KFAv ; ...
      KFBu KFBv ; ...
      KFCu KFCv ] ;

% explicit mixed derivative of the dh/dt misfit term, which depends on B through the thickness h=h(B). It adds to the B row of F^{pq}.
if isB && is_dhdt_meas
    nA=0 ; if isA ; nA=size(KFAu,1) ; end
    iB=nA+(1:size(KFBu,1)) ;
    KFpq(iB,:)=KFpq(iB,:)+JhdotBq(CtrlVar,MUA,F,Meas) ;
end

xi=[KdudA KdudB KdudC ;...
    KdvdA KdvdB KdvdC] ;

if issparse(xi) && nnz(xi)/numel(xi)>0.5  % much faster than sparse times sparse multiplication, unless xi is truly sparse (which it typically will not be)
    xi=full(xi);
end

K=KFpq*xi;

KHess_qp=K+K' ;


end