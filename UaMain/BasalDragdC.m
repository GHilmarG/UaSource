function dbetadC=BasalDragdC(CtrlVar,He,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0)

%% Derivative of the basal drag with respect to the slipperiness C (inverse problem)
%
%   dbetadC=BasalDragdC(CtrlVar,He,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0)
%
% The sliding laws are isotropic, taubx = beta ub and tauby = beta vb (see SlidingLawBeta.m), so that
%
%   dtaubx/dC = dbetadC ub ,   dtauby/dC = dbetadC vb .
%
% Returns dbetadC. The inputs are as for BasalDrag.m (delta and the melange terms are not needed: the melange drag does not depend
% on C).
%
% Replaces the CtrlVar.Inverse.dFuvdClambda=true branch of the earlier BasalDrag.m, where this quantity was returned as the first output
% (taubx). For the Coulomb sliding law, which does not depend on C, an error is raised.
%
% See also: BasalDrag, SlidingLawBeta, SlidingLawCode
%
% (26 Oct 2026)
%%

He(He<0)=0;

if isfield(CtrlVar,"SlidingLawCode")
    Code=CtrlVar.SlidingLawCode;
else
    Code=SlidingLawCode(CtrlVar.SlidingLaw);
end

delta=zeros(size(He),'like',He);   % not needed for dbetadC
[~,~,~,dbetadC]=SlidingLawBeta(Code,CtrlVar,[false false true],He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0);

end
