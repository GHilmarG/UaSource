function W=HessianMultiplyHVP(p,lambda,V,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)  %#ok<INUSL>

%%
%
%   W=HessianMultiplyHVP(p,lambda,V,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)
%
% Hessian-vector products W=H*V at p, used as the HessianMultiplyFcn of fmincon (interior-point algorithm with SubproblemAlgorithm="cg")
% when CtrlVar.Inverse.HessianCalculation="-HVP-". The Hessian is never formed.
%
% The sparse blocks needed (AssembleHessianBlocks.m) are assembled once for each new p, and kept for subsequent products at the same p.
% This requires a forward and an adjoint solve at p (JGH.m), in addition to the one done by fmincon for the cost function and gradient.
% CtrlVar.Inverse.HVP.Type selects the Gauss-Newton ("-GN-", default) or the full ("-full-") Hessian.
%
% lambda (the Lagrange multipliers of nonlinear constraints) is not used, as there are no nonlinear constraints.
%
% see also: AssembleHessianBlocks.m, HessianVectorProduct.m, InversionUsingMatlabOptimisationToolbox3.m
%
%%

persistent pLast Blocks

if isempty(pLast) || numel(pLast)~=numel(p) || any(pLast(:)~=p(:))

    Type="-GN-";
    if isfield(CtrlVar.Inverse,"HVP") && isfield(CtrlVar.Inverse.HVP,"Type")
        Type=CtrlVar.Inverse.HVP.Type;
    end

    CtrlVar.JGH.CalcHessian=false;
    [~,~,~,JGHouts,Fp]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);

    CtrlVarP=CtrlVar ;
    CtrlVarP.Inverse.CholeskyMappingOfCostFunctionAndGradient=false ;
    K=[] ; if isfield(JGHouts,"dFduv") ; K=JGHouts.dFduv ; end   % (16 Nov 2026) Jacobian from the forward solve in JGH.m
    Blocks=AssembleHessianBlocks(CtrlVarP,MUA,Fp,BCs,l,Priors,Meas,BCsAdjoint,JGHouts.Psi_x,JGHouts.Psi_y,Type,K);
    pLast=p;

end

if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient
    % H in the Cholesky-mapped variables: R^{-T} P' H P R^{-1}, see JGH.m
    W=MUA.RG'\(MUA.PRG'*HessianVectorProduct(Blocks,MUA.PRG*(MUA.RG\V)));
else
    W=HessianVectorProduct(Blocks,V);
end

end
