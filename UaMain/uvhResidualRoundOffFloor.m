function rForceFloor=uvhResidualRoundOffFloor(CtrlVar,MUA,F1,Fext0)
%%
% rForceFloor=uvhResidualRoundOffFloor(CtrlVar,MUA,F1,Fext0)
%
% (8 Oct 2026) Estimate of the round-off floor of the uvh cost function rForce (see CalcCostFunctionNRuvh), ie the smallest value
% of rForce that can be reached in floating-point arithmetic for the current state.
%
% With the blockwise normalisation, the h block of the residual is normalised by sh=||Fext0_h||, which is proportional to the time
% step dt (it is the mass flux dt*rho*(|a|+1) integrated over the nodal areas). The round-off error in the h residual, however, is
% given by the representation of the thickness itself, |R_h,i| ~ eps*rho_i*m_i*h_i, and does not decrease with dt. The resulting
% floor of the h contribution to rForce is therefore
%
%    rForceFloor = ( c * eps * ||rho .* m .* |h1| || / sh )^2      (grows as 1/dt^2)
%
% where m are the nodal area weights and c=CtrlVar.uvhResidualRoundOffFloorFactor (default 1, a conservative estimate).
% For small time steps this floor can exceed the desired tolerance (eg 1e-15), and the convergence criterion can then not be
% satisfied. In uvhResidualsCriteria the force tolerance is therefore replaced by max(tolerance,rForceFloor).
% The round-off floor of the momentum blocks is far smaller and is ignored.
%
% Returns 0 if the residual normalisation is not "blockwise".
%%

narginchk(4,4)

rForceFloor=0;
if ~(isfield(CtrlVar,"uvhResidualNormalisation") && CtrlVar.uvhResidualNormalisation=="blockwise")
    return
end

N=MUA.Nnodes;
sh=norm(Fext0(2*N+1:3*N));        % as in CalcCostFunctionNRuvh
if sh<1000*eps
    sh=1;
end

mNode=accumarray(MUA.connectivity(:),repmat(MUA.EleAreas(:)/MUA.nod,MUA.nod,1),[N 1]);   % nodal area weights (positive)

c=1;
if isfield(CtrlVar,"uvhResidualRoundOffFloorFactor")
    c=CtrlVar.uvhResidualRoundOffFloorFactor;
end

rForceFloor=full( (c*eps*norm(F1.rho(:).*mNode.*abs(F1.h(:)))/sh)^2 );

end
