











function [Psi,e]=EnergyDensity(CtrlVar,MUA,F)

narginchk(3,3)


%% 
%
% For the elastic formulation this is strain Energy Density
%
% As explained in, for example:
%
% Vega, B., Yang, J., Tchelepi, H., & Kovscek, A. R. (2020). Investigation of Stress Field and Fracture Development During
% Shale Maturation Using Analog Rock Systems. Transport in Porous Media, 131(2), 503–535.
% https://doi.org/10.1007/s11242-019-01355-2
%
% Griffith's criterion for fracture propagation can be expressed as
%
% $$\partial_l \int \sigma_{ij} \epsilon_{ji} \; dV = \partial_l (G_c A_c) $$
%
% where $l$ is the crack length, $G_c$ is the critical energy release rate, and $A_c$ is the crack surface area.
%
% Franfort and Marigo (1998) proposed another fracture propagation criterion that is able to predict the initiation of a
% fracture as well as fracture path:
%
% $$\int_{V \setminus A_c} \Psi \; dV + \int_{A_c} G_c \; dH^{dim-1} $$
%
%
% The energy $E$ of the system is
%
% $$E = \frac{1}{2} \int_{V \setminus A_c}  \sigma_{ij} \, \epsilon_{ji} \; dV - \int_{\partial V} \mathbf{s}(\mathbf{n}) \cdot \mathbf{d} \; dS + G_c  A_c $$
%
% In the phase-field method, the crack area $A_c$, is approximated as
%
% $$ A_c \approx \frac{1}{2} \int_{V} \left (\frac{\phi^2}{l}  + \| \nabla \phi \|^2 \right ) \; dV $$  
%
% In ice shelves, vertical shear stresses and, in particular, vertical shear strains, are negligible compared to horizontal stresses and
% strains. Therefore
%
% $$ \sigma_{ij} \epsilon_{ij}  \approx \sigma_{xx} \epsilon_{xx} + \sigma_{yy} \epsilon_{yy} + \sigma_{zz} \epsilon_{zz} + 2 \tau_{xy} \epsilon_{xy}$$
% 
%
% Assuming $\mu=1/2$ we have 
% 
% $$\epsilon_{zz}=-(\epsilon_{xx}+ \epsilon_{yy}) $$
% 
%  and
%
% $$ \sigma_{xx} -\sigma_{zz} = 2 \tau_{xx} + \tau_{yy} $$ 
% $$ \sigma_{yy} -\sigma_{zz} = 2 \tau_{yy} + \tau_{xx} $$
%
% leading to
% 
% $$ \sigma_{ij} \epsilon_{ij}  = \sigma_{xx} \epsilon_{xx} + \sigma_{yy} \epsilon_{yy} - \sigma_{zz} (\epsilon_{xx}+\epsilon_{yy}) + 2 \tau_{xy} \epsilon_{xy} $$
% $$  = (\sigma_{xx} -\sigma_{zz} ) \epsilon_{xx} + (\sigma_{yy}-\sigma_{zz}) \epsilon_{yy} +2 \tau_{xy} \epsilon_{xy} $$
% $$  = (2 \tau_{xx} + \tau_{yy} ) \epsilon_{xx} + (\tau_{xx}+2 \tau_{yy}) \epsilon_{yy} + 2 \tau_{xy} \epsilon_{xy} $$
% 
% The work boundary term 
%
% $$ \int_{\partial V} \mathbf{s}(\mathbf{n}) \cdot \mathbf{d} \; dS $$
%
% is calculated in 
% 
%   EdgeWork=EdgeWorkIntegral(CtrlVar,MUA,Displacement,Txy,options)
% 
%
% Using the expression for the vertically integrated traction $\mathbf{s}(\mathbf{n)}$ across the vertical faces of a
% floating ice shelf:
%
% $$ \mathbf{s}(\mathbf{n}) = \frac{1}{4} \rho g (1-\rho/\rho_w) h \, \mathbf{n} $$ 
%
%
%%



%% Calculating internal energy density function 

[dudx,dudy]=calcFEderivativesMUA(F.ub,MUA,CtrlVar) ; 
[dvdx,dvdy]=calcFEderivativesMUA(F.vb,MUA,CtrlVar) ; 
[dudx,dudy,dvdx,dvdy]=ProjectFintOntoNodes(CtrlVar,MUA,dudx,dudy,dvdx,dvdy);

exx=dudx;
eyy=dvdy;
exy=0.5*(dudy+dvdx);
% ezz=-(exx+eyy) ;

e=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));

K=1./(2*F.AGlen);   % Here I use the degraded A (i.e. shear modulus, G) 
txx=2*K.*exx;
tyy=2*K.*eyy;
%tzz=2*K.*ezz;
txy=2*K.*exy; 

Psi=(2*txx+tyy).*exx+(txx+2*tyy).*eyy+2*txy.*txy ; 
Psi(Psi<0)=0; 

return

%% Boundary traction calculations:
% As far as I can see, calculating the work of boundary tractions is not needed
Txy= 0.25*F.rho.*F.g.*(1-F.rho./F.rhow).*F.h  ; 
% Psi=2*A0.^(-1./F.n) .* e.^((F.n+1)./F.n) ; 

Displacement=[F.ub F.vb];
EdgeWork=EdgeWorkIntegral(CtrlVar,MUA,Displacement,Txy,Plots=true); 



end

