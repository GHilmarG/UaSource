function [beta,DbetaU,dbetadh,dbetadC]=SlidingLawBeta(Code,CtrlVar,Need,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0)

%% Drag coefficient beta of the sliding laws, and its derivatives
%
%   [beta,DbetaU,dbetadh,dbetadC]=SlidingLawBeta(Code,CtrlVar,Need,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0)
%
% All sliding laws in Ua are isotropic: the basal drag is
%
%   taubx = beta ub ,   tauby = beta vb ,
%
% with a drag coefficient beta(U,h,C), where U = sqrt(ub^2+vb^2+u0^2) is the (regularised) speed, u0=CtrlVar.SpeedZero. All
% derivatives therefore follow from those of beta (BasalDrag.m):
%
%   dtaubx/dub = beta + DbetaU ub^2 ,   dtaubx/dvb = DbetaU ub vb ,   dtaubx/dh = dbetadh ub ,   dtaubx/dC = dbetadC ub
%   dtauby/dvb = beta + DbetaU vb^2 ,   dtauby/dub = DbetaU ub vb ,   dtauby/dh = dbetadh vb ,   dtauby/dC = dbetadC vb
%
% where DbetaU = (1/U) dbeta/dU.
%
% Inputs:
%
%   Code      : the sliding law (see SlidingLawCode.m)
%   Need      : [NeedU NeedH NeedC], logical: which of DbetaU, dbetadh and dbetadC to calculate (the others are returned as zeros)
%   He, delta : HeavisideApprox and DiracDelta of h-hf (the grounded/floating mask and its derivative with respect to h), at least zero
%   h, H      : ice thickness, and H=S-B (water depth for floating ice)
%   rho, rhow, g : ice and ocean densities, gravitational acceleration
%   ub, vb    : basal velocities
%   C, m      : slipperiness and stress exponent
%   q, muk, V0 : further parameters of some of the sliding laws (exponent of N, Coulomb friction coefficient, threshold speed)
%
% Sliding laws, with Ct=C+C0 (C0=CtrlVar.Czero), t=(U/Ct)^(1/m), b=t/U=Ct^(-1/m) U^(1/m-1), and the effective pressure N of the
% zeroth-order hydrology (N=rho g (h-hf)+Nzero for h>hf, Nzero otherwise):
%
%   1  Weertman  ("-W-")               beta = He b
%   2  Budd      ("-B-","W-N0")        beta = He N^(q/m) b
%   3  Coulomb   ("-C-")               beta = muk N / U
%   4  Tsai      ("minCW-N0")          Weertman or Coulomb, whichever gives the smaller drag
%   5  Cornford  ("rpCW-N0")           beta = muk N He b S^(-1/m),     S = (muk N)^m + He^m U/Ct
%                                      ie |tau| = tauC tauW / (tauC^m + tauW^m)^(1/m), with tauC=muk N and tauW=He Ct^(-1/m) U^(1/m)
%   6  Umbi      ("rCW-N0")            beta = muk N He b / E,          E = muk N + He t
%                                      ie 1/|tau| = 1/tauC + 1/tauW
%   7  Joughin   ("rCW-v0")            beta = He b (U+V0)^(-1/m)
%
% The derivatives with respect to h follow from dHe/dh=delta and dN/dh; for Cornford and Umbi they are written in a form without
% negative powers of He. Coulomb does not depend on C, so dbetadC is not defined for it (an error is raised if requested); for Tsai,
% dbetadC is zero where the Coulomb law applies.
%
% These expressions replace (26 Oct 2026) the local functions Weertman, Budd, Coulomb and N0 of the earlier BasalDrag.m, and the generated rpCWN0.m, rCWN0.m and rCWV0.m.
%
% Written so that it can be used by MATLAB Coder: element-wise operations only (the same code is used for arrays in the vectorised
% assembly and for scalars in the compiled kernel), and all outputs are defined on all paths.
%
% See also: BasalDrag, BasalDragdC, SlidingLawCode
%
% (26 Oct 2026)
%%

NeedU=Need(1); Needh=Need(2); NeedC=Need(3);

u0=CtrlVar.SpeedZero;
U=sqrt(ub.*ub+vb.*vb+u0*u0);
Ct=C+CtrlVar.Czero;

Z=zeros(size(U),'like',U);
DbetaU=Z; dbetadh=Z; dbetadC=Z;

% effective pressure (zeroth-order hydrology), for the laws that use it
if Code>=2 && Code<=6
    [N,dNdh]=EffectivePressureN0(CtrlVar,h,H,rho,rhow,g);
else
    N=Z; dNdh=Z;
end

% Weertman-type factor b = Ct^(-1/m) U^(1/m-1), for all laws except Coulomb (one power)
if Code~=3
    t=(U./Ct).^(1./m);
    b=t./U;
else
    t=Z; b=Z;
end

switch Code

    case 1   % Weertman

        beta=He.*b;
        if NeedU ; DbetaU=(1./m-1).*beta./(U.*U) ; end
        if Needh ; dbetadh=delta.*b ; end
        if NeedC ; dbetadC=-beta./(m.*Ct) ; end

    case 2   % Budd

        qm=q./m;
        Nqm=N.^qm;
        beta=He.*Nqm.*b;
        if NeedU ; DbetaU=(1./m-1).*beta./(U.*U) ; end
        if Needh ; dbetadh=delta.*Nqm.*b+He.*qm.*N.^(qm-1).*dNdh.*b ; end
        if NeedC ; dbetadC=-beta./(m.*Ct) ; end

    case 3   % Coulomb

        beta=muk.*N./U;
        if NeedU ; DbetaU=-beta./(U.*U) ; end
        if Needh ; dbetadh=muk.*dNdh./U ; end
        if NeedC
            error('SlidingLawBeta:CoulombC','The Coulomb sliding law does not depend on C, inversion for C is not possible with it.')
        end

    case 4   % Tsai: Weertman or Coulomb, whichever gives the smaller drag

        betaW=He.*b;
        betaC=muk.*N./U;
        s2=ub.*ub+vb.*vb;
        isC=betaC.*betaC.*s2 < betaW.*betaW.*s2;   % |tau_Coulomb|^2 < |tau_Weertman|^2
        beta=betaW.*~isC+betaC.*isC;   % (element-wise selection, exact for finite values; logical indexing of the right-hand side is not used, for MATLAB Coder)
        if NeedU
            DbetaU=(1./m-1).*betaW./(U.*U).*~isC-betaC./(U.*U).*isC;
        end
        if Needh
            dbetadh=delta.*b.*~isC+muk.*dNdh./U.*isC;
        end
        if NeedC
            dbetadC=-betaW./(m.*Ct).*~isC;
        end

    case 5   % Cornford (regularised Coulomb)

        muN=muk.*N;
        muNm=muN.^m;
        Hem=He.^m;
        S=muNm+Hem.*U./Ct;
        Sm=S.^(-1./m);
        beta=muN.*He.*b.*Sm;
        if NeedU ; DbetaU=beta.*((1./m-1)./(U.*U)-Hem./(m.*Ct.*S.*U)) ; end
        if Needh ; dbetadh=b.*Sm./S.*(muNm.*muN.*delta+muk.*dNdh.*Hem.*He.*U./Ct) ; end
        if NeedC ; dbetadC=-beta.*muNm./(m.*Ct.*S) ; end

    case 6   % Umbi (reciprocal Coulomb-Weertman)

        muN=muk.*N;
        E=muN+He.*t;
        beta=muN.*He.*b./E;
        if NeedU ; DbetaU=beta.*((1./m-1)./(U.*U)-He.*b./(m.*E.*U)) ; end
        if Needh ; dbetadh=muk.*b.*(dNdh.*He.*He.*t+delta.*muN.*N)./(E.*E) ; end
        if NeedC ; dbetadC=-beta.*muN./(m.*Ct.*E) ; end

    case 7   % Joughin

        w=(U+V0).^(-1./m);
        beta=He.*b.*w;
        if NeedU ; DbetaU=beta.*((1./m-1)./(U.*U)-1./(m.*(U+V0).*U)) ; end
        if Needh ; dbetadh=delta.*b.*w ; end
        if NeedC ; dbetadC=-beta./(m.*Ct) ; end

    otherwise

        error('SlidingLawBeta:UnknownCode','Unknown sliding-law code %d.',Code)

end

end

%% local functions

function [N,dNdh]=EffectivePressureN0(CtrlVar,h,H,rho,rhow,g)

% Effective pressure of the zeroth-order hydrology: N = rho g (h-hf) + Nzero where h>hf, and Nzero where h<=hf (floating),
% with the flotation thickness hf=rhow H/rho (zero where H<0).

hf=rhow.*H./rho;
hf(hf<eps)=0;
Dh=h-hf;
I=Dh<eps;
Dh(I)=0;

N=rho.*g.*Dh+CtrlVar.Nzero;
dNdh=rho.*g+zeros(size(Dh),'like',Dh);
dNdh(I)=0;

end
