

function rhoE=DensityEffectivePFF(CtrlVar,rhoi0,rhow,phi)

narginchk(4,4)


gphi=DegradationFunction(CtrlVar,phi) ;

rhoE=gphi.*rhoi0+(1-gphi).*rhow ;


end
