function [ah,rh,mLump]=ThicknessConstraintReactions(CtrlVar,MUA,F,BCs,l)

%% Discrete reactions of the thickness constraints, and the equivalent local mass-balance rates
%
%   [ah,rh,mLump]=ThicknessConstraintReactions(CtrlVar,MUA,F,BCs,l)
%
% (15 Oct 2026) Moved here from ActiveSetUpdate.m, so that it can also be used to provide the reactions as an output (F.ah).
%
% rh    : discrete nodal reactions of the thickness constraints, rh=Lh'*lambda_h. These are the terms that the constraints add to the
%         discrete mass-conservation equations (units of mass, as these equations are multiplied by rho and dt). Non-zero only at
%         constrained nodes.
% ah    : the equivalent local mass-balance rate (units of distance/time, eg m/yr), ah=-rh./(rho mLump dt). Positive where the constraint
%         adds ice, negative where it removes ice.
% mLump : HRZ-lumped nodal weights (areas). These are positive also for higher-order elements, and equal to the row sums of the
%         consistent mass matrix for linear elements.
%
% All thickness constraints are included, ie the active-set constraints (h>=ThickMin), the level-set constraints (downstream of calving
% fronts), and user-defined thickness boundary conditions (BCs.hFixedNode).

MLC=BCs2MLC(CtrlVar,MUA,BCs);

if isempty(MLC.hL) || isempty(l.h)
    rh=zeros(MUA.Nnodes,1) ;
else
    rh=full(MLC.hL'*l.h) ;                          % discrete nodal reactions (non-zero only at constrained nodes)
end

MeDiag=zeros(MUA.Nele,MUA.nod) ; EleArea=zeros(MUA.Nele,1) ;
for Iint=1:MUA.nip                               % HRZ-lumped nodal weights (positive also for higher-order elements)
    funI=shape_fun(Iint,2,MUA.nod,MUA.points) ;
    detJwI=MUA.DetJ(:,Iint)*MUA.weights(Iint) ;
    MeDiag=MeDiag+detJwI.*(funI.^2).' ;
    EleArea=EleArea+detJwI ;
end
mEle=MeDiag.*(EleArea./sum(MeDiag,2)) ;
mLump=accumarray(MUA.connectivity(:),mEle(:),[MUA.Nnodes 1]) ;
ah=-rh./(F.rho.*mLump.*F.dt) ;                  % nodal array, units of mass balance (distance/time)

end
