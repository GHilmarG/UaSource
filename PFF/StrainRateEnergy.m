function [Psi,e,eInt]=StrainRateEnergy(CtrlVar,MUA,F,A0)



[dudx,dudy,xint,yint]=calcFEderivativesMUA(F.ub,MUA,CtrlVar) ; 
[dvdx,dvdy]=calcFEderivativesMUA(F.vb,MUA,CtrlVar) ; 

exx=dudx;
eyy=dvdy;
exy=0.5*(dudy+dvdx);

eInt=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));

e=ProjectFintOntoNodes(CtrlVar,MUA,eInt);
e(e<0)=0; 


Psi=2*A0.^(-1./F.n) .* e.^((F.n+1)./F.n) ; 



end

