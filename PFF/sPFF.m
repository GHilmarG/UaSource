



function [sE,hE]=sPFF(CtrlVar,S,b0,rhoi0,rhow,phi)


% Calculates upper ice surface elevation, s, from flotation, given lower surface and densities.
%
%
% Does not conserve ice thickness!

narginchk(6,6)


% rhoE=DensityEffectivePFF(CtrlVar,rhoi,rhow,phi) ;



gphi=DegradationFunction(CtrlVar,phi) ;

% CtrlVar.PhaseFieldFracture.RiftsAre="-thin ice above inviscid water-"; 

switch CtrlVar.PhaseFieldFracture.RiftsAre

    case "-thin ice above inviscid water-"

        h0=(S-b0).*rhow./rhoi0; 
        hE=gphi.*h0+(1-gphi)*CtrlVar.ThickMin ;
        sE=hE+b0;

    case "-viscous water columns-"

        % flotation:  (s-b) rhoE = (S-b) rhow
        hE=(S-b0).*rhow./rhoE; 
        sE=hE+b0 ;
   

    otherwise

        error("case not found")

end

