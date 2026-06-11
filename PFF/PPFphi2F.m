


function  F=PPFphi2F(CtrlVar,MUA,F)

%%  Maps phase field variable (F.phi) to A, n, rho, s, h and GF  
%
%
% $\phi$ is the phase-field variable.
%
% $\phi=0$ for undamaged material 
% $\phi=1$ for (fully) damaged material 
%
%
% $A$ and $\rho$
%
% $$A=A_0/g_{\phi}^n $$
%
% $$\rho=g_{\phi} \rho_0 + (1-g_{\phi} ) \rho_w $$
%
%
% Then there are two options for the ice thickness $h$
%
%
% Thin ice above inviscid water:
%
% $$h=g_{\phi} h_0 $$ 
%
%
% or as viscous water columns:
%
% $$ h=(S-b_0) \rho_w/\rho $$
%
% where $g_{\phi}$ is the degradation function.
%
% $$ g_{\phi}=(1-k) (1-\phi)^2 + k $$
%
% where $k$ is a regularization parameter.
%
%
% use:
%
% F.rho0
% F.AGlen0
% F.b0
% F.h0
%
%
%%

% Shear modulus G
% Young modulus E
% Poisson's ratio mu
%
%  epsilon = tau/(2 G)
%
% G= E/(2*(1+mu)) = E/(2*(1+1/2)) = E/3

nargoutchk(1,1)


if CtrlVar.PhaseFieldFracture.Phase=="-elastic-"

    % Solving the elastic problem to determine energy release and the resulting phase field 
    % E=9 GPa = 9e6 kPa
    % This is a work in progress, and at the moment the value of the elastic shear modulus is hardwired into the code at this
    % place,
    E=9.5e6 ;
    mu=0.5 ; 
    G= E/(2*(1+mu)) ;
    A0=G+zeros(MUA.Nnodes,1);
    n=1+zeros(MUA.Nnodes,1);


else

    A0=F.AGlen0;
    n=F.n;

end

phi=F.phi;

rhoi0=F.rho0;
h0=F.h0;
b0=F.b0;
rhow=F.rhow ;

S=F.S;
B=F.B;

gphi=DegradationFunction(CtrlVar,phi) ;



switch CtrlVar.PhaseFieldFracture.RiftsAre

    case "-thin ice above inviscid water-"

        % Reduced thickness, where afloat.
        % In the limit of gphi=0 (full damage) the ice thickness is ThickMin where afloat.
        % h0=(S-b0).*rhow./rhoi0;

        % Here ONLY the ice thickness is changed!
        hE=gphi.*h0+(1-gphi)*CtrlVar.ThickMin ;
        rhoE=F.rho;
        AE=F.AGlen; 

    case "-viscous water columns-"

        % flotation:  (s-b) rhoE = (S-b) rhow
        hE=(S-b0).*rhow./rhoiE;  % for fully damaged ice, rhoE=rhow and h0=S-b0
        rhoE=gphi.*rhoi0+(1-gphi).*rhow ;
        AE=A0./ (gphi.^n)  ;


    otherwise

        error("case not found")

end

% The idea is that ice thickness over grounded areas is never changed.

hE=h0.*F.GF.node + (1-F.GF.node).*hE; 

[bE,sE,hE,GF]=Calc_bs_From_hBS(CtrlVar,MUA,hE,S,B,rhoE,rhow) ; 


%% Here I do all of the modifications to F
F.AGlen=AE;
F.n=n; 
F.rho=rhoE;
F.s=sE;
F.b=bE;
F.h=hE;
F.GF=GF;
%%

end