


function  [AE,rhoE,sE,bE,hE,GF]=PPFphi2F(CtrlVar,MUA,F)

% use:
%
% F.rho0
% F.AGlen0
% F.b0
% F.h0
%
%

phi=F.phi;
A0=F.AGlen0;
rhoi0=F.rho0;
h0=F.h0;
b0=F.b0;
rhow=F.rhow ; 
n=F.n; 
S=F.S;
B=F.B;

gphi=DegradationFunction(CtrlVar,phi) ;

rhoE=gphi.*rhoi0+(1-gphi).*rhow ;
AE=A0./ (gphi.^n)  ;

switch CtrlVar.PhaseFieldFracture.RiftsAre

    case "-thin ice above inviscid water-"

        % Reduced thickness, where afloat.
        % In the limit of gphi=0 (full damage) the ice thickness is ThickMin where afloat.
        % h0=(S-b0).*rhow./rhoi0;
        hE=gphi.*h0+(1-gphi)*CtrlVar.ThickMin ;
     

    case "-viscous water columns-"

        % flotation:  (s-b) rhoE = (S-b) rhow
        hE=(S-b0).*rhow./rhoiE;  % for fully damaged ice, rhoE=rhow and h0=S-b0 
     
     

    otherwise

        error("case not found")

end

% The idea is that ice thickness over grounded areas is never changed.

hE=h0.*F.GF.node + (1-F.GF.node).*hE; 

[bE,sE,hE,GF]=Calc_bs_From_hBS(CtrlVar,MUA,hE,S,B,rhoE,rhow) ; 

end