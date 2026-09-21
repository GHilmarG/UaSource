

function [eta,E,eEff,detadA,d2etadAdA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy)

narginchk(6,6)
nargoutchk(1,5)


if nargout==1
    % if the number of output arguments is only 2, then clearly there is not need to calculate any derivatives

    CtrlVar.EffectiveViscosity.CalculateDerivatives=false;

else

    % has on call, calculating derivatives been set? If not then set it to true if number of output arguments sufficient
    if ~(isfield(CtrlVar,"EffectiveViscosity") && isfield(CtrlVar.EffectiveViscosity,"CalculateDerivatives"))

        CtrlVar.EffectiveViscosity.CalculateDerivatives=true;
    end
end

eEff=sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2);

ninv=1./n;

%if isscalar(n) || all(n(:)==n(1))  ( the all(n(:)) is too expensive
% scalar exponent: much faster path
% s=real((eEff./A).^(1/n(1)))./eEff;
%else
eta2=real((eEff./A).^ninv)./eEff;  % this is pretty much just 2 eta
%end

eta=0.5*eta2+CtrlVar.etaZero ;

if CtrlVar.EffectiveViscosity.CalculateDerivatives
    E=(1-n)./(4*n).*eta2./eEff.^2;
else
    E=[];
end

if nargout>3
    detadA   = -eta2./(2*n.*A);
    if nargout>4
        d2etadAdA= 0.5*(ninv.^2+ninv).*eta2./A.^2 ;
    end
end


end
