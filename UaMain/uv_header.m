%% Relates to the solution for the velocity components u and v of the system:
%
% 
% 
% $$
% F^x_i= \left \langle  h \eta \, ( 4 \partial_x u + 2 \partial_y v) \vert  \, \partial_x \phi_i \right \rangle
% + \langle   h \eta \, (\partial_y u + \partial_x v)  \vert  \partial_y \phi_i \rangle
% + \langle \mathcal{G} \beta^2\, u , \phi_i \rangle 
%  - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \Big\vert \partial_x \phi_i \right \rangle
% + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_x B \vert  \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0 
% $$
%
% $$
% F^y_i= \langle  h \eta \, ( 4 \partial_y v + 2 \partial_x u) \vert \partial_y \phi_i \rangle
% +\langle   h \eta \, (\partial_x v + \partial_y u)  \vert \, \partial_x \phi_i \rangle 
% + \langle \mathcal{G} \, \beta^2 \, v \vert  \phi_i \rangle 
%   - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) \Big|   \, \partial_y \phi_i \right \rangle
% +  \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_y B \vert \phi_i \rangle=0
% $$
%
%
%
% Here we use
% 
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
% $\mathcal{G}$ is the floating mask, 1 if grounded, 0 if afloat.
%
% $$\mathcal{G}=\mathcal{H}(h-h_f) $$
%
% where $\mathcal{H}$ is the Heaviside step function and
%
%
% $$h_f=(S-B) \rho_o/\rho $$
%
% where:
% 
% $h$ is the ice thickness
%
% $\rho$ the ice density
%
% $\rho_o$ the ocean density 
%
% $B$ the bedrock
%
% $s$ the upper glacier surface
%
% $b$ the lower glacier surface
%
% $S$ the ocean surface
%
% $$\alpha$$ the slope of the vertical axis of the coordinate system with respect to gravity
%
% $u$ the $x$ velocity component
%
% $v$ the $y$ velocity component
%
%
% The effective viscosity is: 
%
% $$
% \eta= \frac{1}{2} A^{-1/n} \, \left ((\partial_x u)^2 + (\partial_y v)^2 + \partial_x u \,\partial_y v + (\partial_x v + \partial_y u)^2/4+\epsilon_0^2 \right)^{(1-n)/2n} +\eta_0
% $$
%
%
% The effective viscosity is therefore a function of the velocity components and the rheological parameters $A$ and $n$.
%
% The function
% 
%   EffectiveViscositySSTREAM.m
%
% returns the effective viscosity, eta, as well as some derivatives with respect to $A$.
%
% In the particular case of Weertman sliding law $$\beta^2$$ is given by:
%
% $$
% \beta^2=(C+C_0)^{-1/m} \; \left (u_b^2+v_b^2+u_0^2 \right)^{(1-m)/2m} 
% $$
%
% $$\beta^2$$ is therefore a function of the velocity components, and the basal sliding law parameter $C$. For more general sliding
% laws $\beta^2$ may depend on some other parameters as well.  
%
% Often the basal drag term is written as
%
%
% $$t_{bx} =\mathcal{G} \beta^2\, u $$
%
% $$t_{by} =\mathcal{G} \beta^2\, v $$
%
%
% where $t_{bx}$ and $t_{by}$ are the basal traction components.
%
% The function
% 
%   BasalDrag.m
%
% returns $t_{bx}$ and $t_{by}$ as well as various derivatives with respect to $u$, $v$,  $h$ and $C$
%
%%