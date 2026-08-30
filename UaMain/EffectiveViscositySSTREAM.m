function [eta,E,eEff,detadA,d2etadAdA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy)

narginchk(6,6)
nargoutchk(1,5)

%%
% Refactored version of EffectiveViscositySSTREAM.m.  Mathematically
% identical for A>0.
%
% All five outputs are built from a single vector-exponent power.  Writing
%
%    s = A^(-1/n) eEff^(1/n-1) = (eEff/A)^(1/n) / eEff
%
% and using (1-n)/n = 1/n-1, (1-3n)/n = 1/n-3, -1/n-1 = -(1/n)-1 and
% -1/n-2 = -(1/n)-2, the outputs are
%
%    eta       = 0.5 s + etaZero
%    E         = (1-n)/(4n) s / eEff^2
%    detadA    = -s / (2 n A)
%    d2etadAdA = 0.5 (1/n^2 + 1/n) s / A^2
%
% The original evaluates four vector-exponent powers for nargout<=3 and six
% for nargout==5.  A power with an array exponent goes through log and exp
% per element and is by far the most expensive operation here, so this is
% close to a fourfold reduction in the dominant cost.
%
% Two further points:
%
%  1) If n is spatially uniform - which it is in most runs - the exponent
%     is a scalar and MATLAB takes a considerably faster path.  This is
%     detected here rather than assumed.
%
%  2) The original applies real() four to six times.  eEff is real by
%     construction, since
%
%        exx^2 + eyy^2 + exx eyy = (exx + eyy/2)^2 + (3/4) eyy^2 >= 0
%
%     so the only way a complex value can arise is A<0, which is
%     unphysical and which the callers exclude by clamping A to
%     CtrlVar.AGlenmin.  real() is therefore applied once, to s, rather
%     than to every output.  Since real(z*r) = real(z)*r for real r, this
%     is equivalent for A>0.
%
%     Note for A<0 the two versions are NOT equivalent, because
%     (eEff/A)^(1/n) and A^(-1/n)*eEff^(1/n-1) lie on different branches.
%     Callers clamp A, so this should not arise, but it is a genuine
%     difference in behaviour.
%
%%

eEff=sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2);

ninv=1./n;

if isscalar(n) || all(n(:)==n(1))
    % scalar exponent: much faster path
    s=real((eEff./A).^(1/n(1)))./eEff;
else
    s=real((eEff./A).^ninv)./eEff;
end

eta=0.5*s+CtrlVar.etaZero ;
E=(1-n)./(4*n).*s./eEff.^2;

if nargout>3
    detadA   = -s./(2*n.*A);
    if nargout>4
        d2etadAdA= 0.5*(ninv.^2+ninv).*s./A.^2 ;
    end
end

end
