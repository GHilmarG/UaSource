




function [AE,rhoE]=ArhoPFF(CtrlVar,phi,rhoi0,rhow,A0,n)

narginchk(6,6)


% calculates AGlen effective and the effective density as a function of phi

gphi=DegradationFunction(CtrlVar,phi) ;
rhoE=DensityEffectivePFF(CtrlVar,rhoi0,rhow,phi) ;

AE=A0./ (gphi.^n)  ;


end



