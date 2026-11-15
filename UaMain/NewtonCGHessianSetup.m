function Hv=NewtonCGHessianSetup(x,Fx,JGHouts,CtrlVar,MUA,BCs,l,Priors,Meas,BCsAdjoint)

%%
%
%   Hv=NewtonCGHessianSetup(x,Fx,JGHouts,CtrlVar,MUA,BCs,l,Priors,Meas,BCsAdjoint)
%
% Returns a function handle, Hv(V), for Hessian-vector products with respect to the optimisation variables x, as used by
% UaOptimisationNewtonCG.m. Fx is the state at x, and JGHouts the additional outputs of JGH.m at x (adjoint variables and gradients).
%
% The blocks are assembled in the physical parameters p (AssembleHessianBlocks.m). With the box transform, p=T(x) element by element,
% and with D=diag(dp/dx),
%
%   Gauss-Newton (CtrlVar.Inverse.HVP.Type="-GN-", default) :  H_x v = D H_p (D v)
%   full         (CtrlVar.Inverse.HVP.Type="-full-")          :  H_x v = D H_p (D v) + (g_p .* d^2p/dx^2) .* v
%
% where g_p is the gradient with respect to p. The Gauss-Newton form omits the term due to the curvature of the box transform (which
% is multiplied by the gradient), and is positive semi-definite.
%
% see also: UaOptimisationNewtonCG.m, AssembleHessianBlocks.m, HessianVectorProduct.m, BoxTransform.m
%
%%

Type="-GN-";
if isfield(CtrlVar.Inverse,"HVP") && isfield(CtrlVar.Inverse.HVP,"Type")
    Type=CtrlVar.Inverse.HVP.Type;
end

CtrlVarP=CtrlVar ;
CtrlVarP.Inverse.BoxTransform=false ;
CtrlVarP.Inverse.CholeskyMappingOfCostFunctionAndGradient=false ;
Blocks=AssembleHessianBlocks(CtrlVarP,MUA,Fx,BCs,l,Priors,Meas,BCsAdjoint,JGHouts.Psi_x,JGHouts.Psi_y,Type);

if CtrlVar.Inverse.BoxTransform

    D=BoxTransform("jacobian",x,MUA.BoxTransform);
    if Type=="-full-"
        gp=JGHouts.dRdp(:)+JGHouts.dIdp(:);
        E=gp.*BoxTransform("curvature",x,MUA.BoxTransform);
    else
        E=zeros(size(D));
    end
    Hv=@(V) D.*HessianVectorProduct(Blocks,D.*V)+E.*V ;

else

    Hv=@(V) HessianVectorProduct(Blocks,V) ;

end

end
