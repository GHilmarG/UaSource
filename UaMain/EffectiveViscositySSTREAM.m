function [eta,E,eEff,detadA,d2etadAdA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy)

narginchk(6,6)
nargoutchk(1,5)



eEff=sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2);

ninv=1./n;

%if isscalar(n) || all(n(:)==n(1))  ( the all(n(:)) is too expensive
% scalar exponent: much faster path
% s=real((eEff./A).^(1/n(1)))./eEff;
%else
s=real((eEff./A).^ninv)./eEff;
%end

eta=0.5*s+CtrlVar.etaZero ;

if CtrlVar.EffectiveViscosity.CalculateDerivatives
    E=(1-n)./(4*n).*s./eEff.^2;
else
    E=[];
end

if nargout>3
    detadA   = -s./(2*n.*A);
    if nargout>4
        d2etadAdA= 0.5*(ninv.^2+ninv).*s./A.^2 ;
    end
end


end
