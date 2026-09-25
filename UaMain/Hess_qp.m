






function [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC)

narginchk(13,13)
nargoutchk(1,1)

%% Builds the mixed Hessian term
%
% $$ H^{pq}+H^{qp} = \big( \mathcal{F}^{pq} \big) \xi + \Big[ \big( \mathcal{F}^{pq} \big) \xi \Big]^T $$
%
% The $J^{pq}$ contribution is absent because each term in the cost function is an explicit function of either p or
% q, but not of both.
%
% Note: that statement relies on the measurements being velocities only. If dh/dt measurements are included then
% $J_{\dot{h}}$ depends on the velocities AND on the thickness, and hence on B, so $J^{pq} \neq 0$ for B. That term
% is not implemented.
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


if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)

    error("Hess_qp:CaseNotImplemented","Hessian mixed qp terms not fully implemented for use the dh/dt meas. ")

end


[KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

[KFBu,KFBv]=FBuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

[KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;


KFpq=[KFAu KFAv ; ...
      KFBu KFBv ; ...
      KFCu KFCv ] ;

xi=[KdudA KdudB KdudC ;...
    KdvdA KdvdB KdvdC] ;

if issparse(xi) && nnz(xi)/numel(xi)>0.5  % much faster than sparse times sparse multiplication, unless xi is truly sparse (which it typically will not be)
    xi=full(xi);
end

K=KFpq*xi;

KHess_qp=K+K' ;


end