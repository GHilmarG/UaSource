function [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh] = ...
    BasalDrag(CtrlVar,Derivatives,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0,uo,vo,Co,mo,ua,va,Ca,ma)

%% Basal drag and its derivatives with respect to ub, vb and h
%
%   [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh] = ...
%       BasalDrag(CtrlVar,Derivatives,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0,uo,vo,Co,mo,ua,va,Ca,ma)
%
% Derivatives : which derivatives to calculate
%                 0 : none (only taubx and tauby)
%                 1 : with respect to ub and vb
%                 2 : with respect to ub, vb and h
%               Derivatives that are not calculated are returned as zeros. All outputs are always defined.
%
% Outputs:
%
%   taubx, tauby       : basal drag
%   dtaubxdu, dtaubxdv : derivatives of taubx with respect to ub and vb
%   dtaubydu, dtaubydv : derivatives of tauby with respect to ub and vb
%   dtaubxdh, dtaubydh : derivatives of taubx and tauby with respect to h
%
%   [ Delta taubx ]  =  [ dtaubxdu   dtaubxdv   dtaubxdh ]  [ Delta ub ]
%   [ Delta tauby ]     [ dtaubydu   dtaubydv   dtaubydh ]  [ Delta vb ]
%                                                           [ Delta h  ]
%
% Inputs:
%
%   CtrlVar   : uses SlidingLaw (or SlidingLawCode, if present), Czero, SpeedZero, Nzero and IncludeMelangeModelPhysics
%   He, delta : HeavisideApprox(kH,h-hf,Hh0) and DiracDelta(kH,h-hf,Hh0), hf=rhow H/rho
%   h, H      : ice thickness, and H=S-B
%   rho, rhow, g : ice and ocean densities, gravitational acceleration
%   ub, vb    : basal velocities
%   C, m      : slipperiness and stress exponent
%   q, muk, V0 : parameters of some of the sliding laws (only used by those)
%   uo, vo, Co, mo, ua, va, Ca, ma : ocean and air velocities and drag parameters of the melange terms (only used if
%               CtrlVar.IncludeMelangeModelPhysics)
%
% The sliding laws are isotropic, taubx = beta ub, tauby = beta vb; the drag coefficient beta and its derivatives are calculated by
% SlidingLawBeta.m for all sliding laws, and the derivatives of the drag follow from those in the same way for all of them.
% If CtrlVar.IncludeMelangeModelPhysics, the drag of the ocean and of the air on the melange, (1-He) beta_o (ub-uo), (1-He) beta_a
% (ub-ua), with Weertman-type drag coefficients, is added.
%
% The derivative of the drag with respect to C (inverse problem) is calculated by BasalDragdC.m.
%
% Rewritten 26 Oct 2026. Differences from the earlier version of BasalDrag.m:
%   - which derivatives are calculated is given explicitly (Derivatives), not by the number of output arguments or by
%     CtrlVar.BasalDrag.CalculateDerivatives
%   - the derivative with respect to C is a separate function (BasalDragdC.m), instead of being returned in taubx if
%     CtrlVar.Inverse.dFuvdClambda=true
%   - the names and aliases of the sliding laws are defined in SlidingLawCode.m only ("-W-", "-B-", "-C-"; "W", "B", "C" no longer)
%   - MUA and B are not inputs, and the melange drag is not returned separately
%
% Written so that it can be used by MATLAB Coder: element-wise operations only (arrays in the vectorised assembly, scalars in the
% compiled kernel), all outputs defined on all paths, no string arrays.
%
% See also: SlidingLawBeta, BasalDragdC, SlidingLawCode
%
% (26 Oct 2026)
%%

% rounding and the interpolation to integration points can make these slightly negative
He(He<0)=0;
delta(delta<0)=0;

if isfield(CtrlVar,"SlidingLawCode")
    Code=CtrlVar.SlidingLawCode;     % eg set once per assembly (compiled kernel)
else
    Code=SlidingLawCode(CtrlVar.SlidingLaw);
end

if coder.target('MATLAB')
    if any(C+CtrlVar.Czero<0,'all')
        fprintf("BasalDrag: C+CtrlVar.Czero is negative, min(C)=%g. A possible solution might be to increase CtrlVar.Czero. \n",min(C,[],'all'))
    end
end

Need=[Derivatives>=1, Derivatives>=2, false];
[beta,DbetaU,dbetadh]=SlidingLawBeta(Code,CtrlVar,Need,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0);

taubx=beta.*ub;
tauby=beta.*vb;

Z=zeros(size(taubx),'like',taubx);
dtaubxdu=Z; dtaubxdv=Z; dtaubydu=Z; dtaubydv=Z; dtaubxdh=Z; dtaubydh=Z;

if Derivatives>=1
    dtaubxdu=beta+DbetaU.*ub.*ub;
    dtaubydv=beta+DbetaU.*vb.*vb;
    dtaubxdv=DbetaU.*ub.*vb;
    dtaubydu=dtaubxdv;
end

if Derivatives>=2
    dtaubxdh=dbetadh.*ub;
    dtaubydh=dbetadh.*vb;
end

if CtrlVar.IncludeMelangeModelPhysics
    % drag of the ocean (uo,vo) and of the air (ua,va) on the melange
    [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh]=AddRelativeDrag(CtrlVar,Derivatives,He,delta,ub,vb,uo,vo,Co,mo, ...
        taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh);
    [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh]=AddRelativeDrag(CtrlVar,Derivatives,He,delta,ub,vb,ua,va,Ca,ma, ...
        taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh);
end

end

%% local functions

function [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh]=AddRelativeDrag(CtrlVar,Derivatives,He,delta,ub,vb,ur,vr,Cr,mr, ...
    taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh)

% Adds the drag (1-He) beta_r (ub-ur, vb-vr) of a medium moving with velocity (ur,vr) (ocean or air acting on the melange), with
% the Weertman-type drag coefficient beta_r = (Cr+C0)^(-1/mr) Ur^(1/mr-1), Ur = sqrt((ub-ur)^2+(vb-vr)^2+u0^2).

U=ub-ur;
V=vb-vr;
Ur=sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2);
Crt=(Cr+CtrlVar.Czero).^(-1./mr);
br=Crt.*Ur.^(1./mr-1);
G=1-He;

taubx=taubx+G.*br.*U;
tauby=tauby+G.*br.*V;

if Derivatives>=1
    Dbr=(1./mr-1).*Crt.*Ur.^(1./mr-3);
    dtaubxdu=dtaubxdu+G.*(br+Dbr.*U.*U);
    dtaubydv=dtaubydv+G.*(br+Dbr.*V.*V);
    dtaubxdv=dtaubxdv+G.*Dbr.*U.*V;
    dtaubydu=dtaubydu+G.*Dbr.*U.*V;
end

if Derivatives>=2
    dtaubxdh=dtaubxdh-delta.*br.*U;
    dtaubydh=dtaubydh-delta.*br.*V;
end

end
