%% Formulation of the momentum equations, using the SSA approximation.
%
% 
% 
% By defining the operators
%
%
% $$ \mathcal{F}^x := \partial_x ( 4 h \eta \partial_x u + 2 h \eta \partial_y v) +\partial_y ( h \eta \, (\partial_x v +
% \partial_y u)) - \mathcal{G} \, \beta^2 \, u -  (\rho g h \partial_x s + \frac{1}{2} h^2 g  \, \partial_x \rho  ) \cos
% \alpha  +  \rho g h \sin \alpha $$
%
% $$
%  \mathcal{F}^y := \partial_y ( 4 h \eta \partial_y v + 2 h \eta \partial_x u) +\partial_x ( h \eta \, (\partial_y u +
%  \partial_x v)) - \mathcal{G} \, \beta^2 \, v
% -  (\rho g h \partial_y s + \frac{1}{2} h^2 g  \, \partial_y \rho  ) \cos \alpha $$
%
% we can write the vertically integrated SSA field equations for momentum in $x$ and $y$ directions as
%
% $$ \mathcal{F}^x = 0 $$
% 
% $$\mathcal{F}^y  = 0 $$
%
%
% The finite-element formulation of this problem is:
%
% $$ F^x_i= \left \langle  h \eta \, ( 4 \partial_x u + 2 \partial_y v) \vert  \, \partial_x \phi_i \right \rangle + \langle
% h \eta \, (\partial_y u + \partial_x v)  \vert  \partial_y \phi_i \rangle + \langle \mathcal{G} \beta^2\, u \vert \phi_i
% \rangle
%  - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \Big\vert \partial_x \phi_i \right \rangle
% + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_x B \vert  \phi_i \rangle  - \langle \rho g \sin(\alpha) \,
% h  | \phi_i \rangle   =0 $$
%
% $$ F^y_i= \langle  h \eta \, ( 4 \partial_y v + 2 \partial_x u) \vert \partial_y \phi_i \rangle +\langle   h \eta \,
% (\partial_x v + \partial_y u)  \vert \, \partial_x \phi_i \rangle + \langle \mathcal{G} \, \beta^2 \, v \vert  \phi_i
% \rangle
%   - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) \Big|   \, \partial_y \phi_i \right \rangle
% +  \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_y B \vert \phi_i \rangle=0 $$
%
% Additionally we have the vertically integrated equation of the conservation of mass which we write as
%
% $$ \mathcal{F}^h = \rho \partial_t h + \partial_x ( \rho h u) + \partial_y (\rho h v) - \rho a =0 $$
%
%
% We solve this using the $\Theta$ method as
%
% $$ \rho \, \frac{1}{\Delta t} (h_1- h_0 ) = (1-\Theta) \left (  \rho\, a(h_0) - \nabla \cdot ( \rho \, \mathbf{v}_0 h_0 )
% \right ) + \Theta \left (  \rho\, a(h_1) - \nabla \cdot ( \rho \, \mathbf{v}_1 h_1 ) \right )  $$
%
% where $\Theta$ is generally set to $\Theta=1/2$
%
%
% Generally the forward problem is non-linear with respect to the variables $u$ and  $v$ and $h$. This is because $\eta$ (the
% effective viscosity), $\beta^2$ (the basal traction coefficient), and $\mathcal{G}$ (the flotation mask) can be functions
% of $u$, $v$ and $h$.
%
% Typically in ice-sheet modeling we have
% 
% $$\eta=\eta(u,v) $$
%
% $$\beta^2=\beta^2(u,v)$$
%
% $$\mathcal{G}=\mathcal{G}(h) $$
%
%
% The $uv$ forward problem is solving the equations
% 
% $$ \mathcal{F}^x = 0 $$
% 
% $$\mathcal{F}^y  = 0 $$
%
% for $u$ and $v$, given $h$, $S$, $B$ and $\rho$, $\rho_o$, $g$, $\alpha$ as well as any material parameters entering the
% expression for the effective viscosity, and the basal sliding law parameters entering the expression for $\beta^2$.
%
% The $uvh$ forward problem is solving the equations for
%
% $$ \mathcal{F}^x = 0 $$
% 
% $$\mathcal{F}^y  = 0 $$
%
% $$\mathcal{F}^h =0 $$
%
% for $u$, $v$ and $h$, given $S$, $B$ and $\rho$, $\rho_o$, $g$, $\alpha$ as well as any material parameters entering the
% expression for the effective viscosity, and the basal sliding law parameters entering the expression for $\beta^2$.
%
%
%
%% Floating relationships, and geometrical closure
%
% Where the ice is afloat:
%
% $$\rho g \, (s-b)  = \rho_o g \, (S-b), $$
%
% $$b=\frac{1}{1-(1-\mathcal{G}) \rho/\rho_o} \left ( \mathcal{G} B + (1-\mathcal{G}) ( S - \rho s/\rho_o ) \right ) $$
%
% $$h= s-\frac{1}{1-(1-\mathcal{G}) \rho/\rho_o} \left ( \mathcal{G} B + (1-\mathcal{G}) ( S - \rho s/\rho_o ) \right ) $$
%
% The geometrical closure condition is
% 
% $$\mathcal{C}:= b - \mathcal{G}(h-h_f) \, (B-b_f) - b_f =0 $$
%
% where
%
% $$b_f:=\frac{\rho s - \rho_o S}{\rho-\rho_o} $$
%
% $$h_f:= \frac{\rho_o}{\rho} (S-B) $$
%
% $h_f$ is the flotation ice thickness. It is the ice thickness exactly at flotation limit, and the ice thickness at the
% grounding line.
%
% $b_f$ is the location of the lower ice surface $b$ when afloat.
%
%
% As a function of $B$, the closure condition $\mathcal{C}$ is a non-linear function as $\mathcal{G}$ depends on $B$, through
% $h_f=h_f(B)$.
%
%
%% Monolithic approach
%
% Using Newton-Raphson we can solve the $uvh$ system using a monolithic approach as
%
% $$ \left( \matrix{ K_{uu}  & K_{uv}  & K_{uh} \cr K_{vu}  & K_{vv}  & K_{vh} \cr K_{hu}  & K_{hv}  & K_{hh}  \cr} \right)
% \left( \matrix{ \Delta u \cr
%                 \Delta v \cr
%                \Delta h} \right)
% = - \left( \matrix{ F^x \cr F^y \cr F^h } \right) $$
%
% This Newton system is then iterated until the norm of the right-hand-side is below a given tolerance.
% 
% The norm needs to be scaled in some sensible way, and here one might use different normalization factors for the momentum
% terms ($F^x$ and $F^y$), and the mass conservation term, ($F^h$).
%
%
% 
%% The staggered approach 
%
% The staggered $uv$ - $h$ solve can be written as:
%
% Choose an initial velocity guess $(u_{1,i},v_{1,i})$ for $i=1$,
%
% For $i=1,2,\ldots $ do
%
% 1) solve the mass conservation equation at $(u_{1,i},v_{1,i})$ for $h_{1,i+1}$ as
%
% $$K_{hh} \, \Delta h^k = - F^h(u_{1,i},v_{1,i},h_{1,i+1}^k ; u_0,v_0, h_0) $$
%
% $$h^{k+1}_{1,i+1}= h^{k}_{1,i+1} + \Delta h^k $$
%
% by iterating over the Newton step number $k$. The $K_{hh}$ block is evaluated at $(u_{1,i},v_{1,i},h^k_{1,i+1})$.
%
% The initial start for the $h$ Newton iteration is $h^0_{1,i+1}=h_{1,i}$.
%
% Once converged, set
%
% $$h_{1,i+1}= h^k_{1,i+1} $$
%
% 2) Solve the momentum equations at $h=h_{1,i+1}$ for $(u_{1,i+1},v_{1,i+1})$
%
% $$ \left( \matrix{ K_{uu}  & K_{uv}   \cr K_{vu}  & K_{vv}   \cr
%  }\right)
% \left( \matrix{ \Delta u^k \cr
%                 \Delta v^k
%                     } \right)
% = - \left( \matrix{
%   F^x(u^k_{1,i+1},v^k_{1,i+1},h_{1,i+1}) \cr F^y(u^k_{1,i+1},v^k_{1,i+1},h_{1,i+1})
%  } \right)
% $$
%
% $$u_{1,i+1}^{k+1} = u^k_{1,i+1}+ \Delta u^k$$
%
% $$v_{1,i+1}^{k+1} = v^k_{1,i+1} + \Delta v^k$$
% 
% by iterating over $k$, Newton iteration number. The $K$ blocks are evaluated at $(u^k_{1,i+1},v^k_{1,i+1},h_{1,i+1})$
%
% The initial start for the $uv$ Newton iteration is $u^0_{1,i+1}=u_{1,i}$, and $v^0_{1,i+1}=v_{1,i}, with the additional
% constraint that the iterate must be feasible. 
%
% Once converged, set
%
% $$u_{1,i+1} = u^k_{1,i+1}$$
%
% $$v_{1,i+1} = v^k_{1,i+1}$$
%
% Unless tolerances (see below) have been met, set $i=i+1$ and go back to step 1
%
%
%% tolerances for the uv-h staggered approach
%
% When solving the staggered system, several convergence criteria are used: For both the individual $uv$ and $h$ solves, the
% convergence criteria are the norms of their respective right-hand sides of the resulting Newton system. And for the
% solution of the combined $uv$-$h$ system the convergence is either based on:
% 
% 1) the changes in the solution vectors, or
%
% 2) by evaluating the norm of the rhs of the full $uvh$ system.
% 
% Both criteria control the full residual, since $F^h$ is proportional to the divergence of the velocity change. Approach 2
% is preferred because it is directly consistent with the monolithic exit criterion and uses the same norm and scaling.
%
% Neither criterion alone guarantees a small error in the solution: when the staggered iteration contracts slowly, small
% changes or small residuals may still correspond to a larger error. 
%
% For approach 2, the $uv$ residual will be within the inner $uv$ tolerance because the momentum solver has been iterated to full convergence already
% for $h_{1,i+1}$. The full residual is then the norm (subject to some scaling) of only
%
% $$F^h(u_{1,i+1},v_{1,i+1},h_{1,i+1}) = \Theta \, \nabla \cdot (\rho \, h_{1,i+1} (\mathbf{v}_{1,i+1} - \mathrm{v}_{1,i}) $$
%
% The velocity change criterion therefore controls the full residual.
%
%% Definitions:
%
%
% The effective viscosity is:
%
% $$ \eta= \frac{1}{2} A^{-1/n} \, \left ((\partial_x u)^2 + (\partial_y v)^2 + \partial_x u \,\partial_y v + (\partial_x v +
% \partial_y u)^2/4+\epsilon_0^2 \right)^{(1-n)/2n} +\eta_0 $$
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
% $$ \beta^2=(C+C_0)^{-1/m} \; \left (u_b^2+v_b^2+u_0^2 \right)^{(1-m)/2m} $$
%
% $$\beta^2$$ is therefore a function of the velocity components, and the basal sliding law parameter $C$ and the stress
% exponent $m$.
% 
% For more general sliding laws $\beta^2$ may depend on some other parameters as well.
%
% The function
%
%   BasalDrag.m
%
% calculates quantities related to the basal drag terms.
%
% Often the basal drag term is written as
%
%
% $$t_{bx} =\mathcal{G} \beta^2\, u $$
%
% $$t_{by} =\mathcal{G} \beta^2\, v $$
%
% and the function
%
%   BasalDrag.m
%
% returns $t_{bx}$ and $t_{by}$ and not $\beta^2$ where $t_{bx}$ and $t_{by}$ are the basal traction components.
%
% 
%
% $\mathcal{G}$ is the flotation mask, 1 if grounded, 0 if afloat.
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
% $h=s-b$ is the ice thickness
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
% $d$ is the submarine ice thickness (always positive), defined as:
% 
% $$d=\mathcal{H}(h_f-h) \, \rho h / \rho_o + \mathcal{H}(h-h_f) \, H^{+} $$
%
% which we can also write as
%
% $$d= (1-\mathcal{G} )\, \frac{\rho}{\rho_o}  h + \mathcal{G} \, H^{+} $$
%
% $$H^{+} = \mathcal{H}(H) \, H $$
%
% $$H=S-B $$
%
%
% Note that because of the identities:
% 
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_x B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_x b $$
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
% we can always replace $B$ in the expressions for $F^x_i$ and $F^y_i$ with $b$.
%
% The function
% 
%   BasalDrag.m
%
% returns $t_{bx}$ and $t_{by}$ as well as various derivatives with respect to $u$, $v$,  $h$ and $C$
%
%%