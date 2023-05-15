function [Eta,E,e,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy)

narginchk(6,6)
nargoutchk(1,4)


e=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));


% What is a reasonable lower limit on viscosity?
%
% 10^15 /1000 /(365.25*24*60*60)  = 31,688 kPa yr
% Viscosity of
% 10^15 /1000 /(365.25*24*60*60)  = 31,688 kPa yr


Eta=real(0.5*A.^(-1./n).*e.^((1-n)./n))+CtrlVar.etaZero ; 

I=Eta>1e8 ;
Eta(I)=1e8; 

%numel(find(Eta<2*eta0))/numel(exx) 


E=real((1-n)./(4*n).*A.^(-1./n).*e.^((1-3*n)./n));



if nargout>3
    dEtadA   = - real(A.^(-1./n-1).*e.^((1-n)./n)./(2*n));
    dEtadA(I)=0 ; 
end

end
