







function dIdA=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)


narginchk(7,7)

%% Calculates the vector quantity:
%
%
% $$ \langle  \delta_{A_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{A_i} F^y \phi_i| \Psi_y \rangle $$
%
% We expand all variables in the same FE basis 
%
% $$ \{ \phi_i \vert i=1 \ldots n \} $$
%
% and arrive at
%
%
% $$ \langle  \partial_A F^x \phi_i | \phi_j \rangle \, \Psi_{x,j} + \langle  \partial_{A} F^y \phi_i| \phi_j \rangle \, \Psi_{y,h}  $$
%
% We do not calculate the matrices 
%
% $$ \langle  \partial_A F^x \phi_i | \phi_j \rangle  $$
% 
% $$ \langle  \partial_A F^y \phi_i| \phi_j \rangle  $$
%
% but assembly the vector directly without the need to store these matrices.
%
%% Formulation of the momentum equations, using the SSA approximation.
%
% 
% 
% By defining the operators
%
%
% $$ 
% \mathcal{F}^x := \partial_x ( 4 h \eta \partial_x u + 2 h \eta \partial_y v) +\partial_y ( h \eta \, (\partial_x v + \partial_y u)) - \mathcal{G} \, \beta^2 \, u 
% -  (\rho g h \partial_x s + \frac{1}{2} h^2 g  \, \partial_x \rho  ) \cos \alpha  +  \rho g \sin \alpha 
% $$
%
% $$
%  \mathcal{F}^y := \partial_y ( 4 h \eta \partial_y v + 2 h \eta \partial_x u) +\partial_x ( h \eta \, (\partial_y u + \partial_x v)) - \mathcal{G} \, \beta^2 \, v  
% -  (\rho g h \partial_y s + \frac{1}{2} h^2 g  \, \partial_x \rho  ) \cos \alpha  
% $$
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
% Additionally we have the vertically integrated equation of the conservation of mass which we write as 
%
% $$ \mathcal{F}^h = \rho \partial_t h + \partial_x ( \rho h u) + \partial_y (\rho h v) - \rho a =0 $$
%
% Generally the forward problem is non-linear with respect to the variables $u$ and  $v$. This is because 
% $\eta$ (the effective viscosity), $\beta^2$ (the basal traction coefficient), and $\mathcal{G}$ (the flotation mask) can be functions of $u$, $v$
% and $h$.
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
% for $u$ and $v$, given $h$, $S$, $B$ and $\rho$, $\rho_o$, $g$, $\alpha$ as well as any material parameters entering
% the expression for the effective viscosity, and the basal sliding law parameters entering the expression for $\beta^2$.
%
% The $uvh$ forward problem is solving the equations for 
%
% $$ \mathcal{F}^x = 0 $$
% 
% $$\mathcal{F}^y  = 0 $$ 
%
% $$\mathcal{F}^h =0 $$
%
% for $u$, $v$ and $h$, given $S$, $B$ and $\rho$, $\rho_o$, $g$, $\alpha$ as well as any material parameters entering
% the expression for the effective viscosity, and the basal sliding law parameters entering the expression for $\beta^2$.
%
%
% The inverse problem is selecting one or more of $A$, $B$ and $C$ given measurements of $u$ and $v$, and $\partial_t h$ while respecting any
% other direct/prior information about $A$, $B$ and $C$. 
%
% For the inversion problem, the forward problem are the two momentum equations, i.e. 
%
% 
% $$ \mathcal{F}^x = 0 $$
% 
% $$\mathcal{F}^y  = 0 $$ 
%
%
% Even when we have measurements of $\partial_t h$, we do NOT include the mass-conservation equation  $$\mathcal{F}^h =0 $$
% as part of the forward model. Rather, we include this misfit term directly in the cost function, $J$, that we minimize. 
%
%
% Solving the inverse problem can be expressed as a minimization problem were we minimize a scalar function $J$ with respect
% to the model parameters $p$ where $p$ can be one or more of $A$, $B$ and $C$. 
%
% $p$: one ore more of $A$, $B$ and $C$ 
%
% $q$: the $x$ and $y$ velocity components $u$ and $v$
%
% The inverse problem is then
%
% $$ \min_p J(q(p),p) $$
%
% where $p$ are the model parameters we want to invert for, generally one or more of $A$, $B$ and $C$, and $q$ are the field
% variables, here these are the two velocity components $u$ and $v$.
%
% The cost function $J$ is generally formulated in the continuous limit as an integral over the computational domain, 
%
% $$ J(q(p),p) = \int \!\!\int \mathcal{J}(q(p),p)  \; dx \, dy $$
%
% To find the minima of $J(q(p),p)$ with respect to $p$ while at the same time fulfilling the equations $\mathcal{F}^x(q) = 0$
% and $\mathcal{F}^y(q) = 0$ of the forward problem we use the method of Lagrange parameters and introduce the extended
% functional 
% 
% $$\mathcal{L} = J(p(q),p) + \langle \Psi_x \vert \mathcal{F}^x \rangle + \langle \Psi_y \vert \mathcal{F}^y \rangle $$
%
%
%% Definitions:
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
% $$\beta^2$$ is therefore a function of the velocity components, and the basal sliding law parameter $C$ and the stress exponent $m$. 
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
% Note the because of the identities: 
% 
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
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
%
%% Derivatives 
%
%
% $$ \delta_A  \eta   = -\frac{1}{2n}  \, A^{-1/n-1} \; e^{(1-n)/(2n)} \; \delta A $$
%
% or
%
% $$ \delta_A  \eta   = \partial_A \eta \; \delta A $$
%
%
% The second derivative is:
%
% $$ \delta^2_{AA}  \eta   = \frac{1}{2} \;  \frac{-1}{n} \; (-1/n-1)   \, A^{-1/n-2} \; e^{(1-n)/(2n)} \; \delta A \, \delta A $$
%
% i.e.
%
% $$ \delta^2_{AA}  \eta   = \frac{1}{2} \;  (1/n^2 + 1/n)   \, A^{-1/n-2} \; e^{(1-n)/(2n)} \; \delta A \, \delta A $$
%
% or
%
% $$ \delta^2_{AA}  \eta   = \partial^2_{AA} \eta  \; \delta A \, \delta A $$
%
% We therefore have
%
% $$ 
% \langle  \delta_{A_i} F^x | \Psi_x \rangle = -\langle h \, (d\eta/dA) \, \phi_i \, ( 4 \partial_x u + 2 \partial_y v) | \partial_x \Psi_x \rangle - \langle h \, (d \eta/dA) \, \phi_i \, (\partial_y u + \partial_x v) | \partial_y \Psi_x \rangle  
% $$
%
% that is
%
% $$ 
% \langle  \delta_{A_i} F^x | \Psi_x \rangle 
% = -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \Psi_x  + (\partial_y u + \partial_x v) \, \partial_y \Psi_x \big ) \, \phi_i \; dx \, dy 
% $$
%
% and
%
% $$ 
% \langle  \delta_{A_i} F^x | \Psi_x \rangle + \langle  \delta_{A_i} F^y | \Psi_y \rangle 
% = -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \Psi_x  + (\partial_y u + \partial_x v) \, \partial_y \Psi_x \big ) \, \phi_i \; dx \, dy 
%   -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \Psi_y  + (\partial_x v + \partial_y u) \, \partial_x \Psi_y \big ) \, \phi_i \; dx \, dy 
% $$
%
%
% which we can write as
%
% $$
% \langle  \delta_{A_i} F^x | \Psi_x \rangle + \langle  \delta_{A_i} F^y | \Psi_y \rangle 
% = -\int \, \partial_A \eta \,  \, h \, 
%   \big  (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \Psi_x  + (\partial_y u + \partial_x v) \, \partial_y \Psi_x 
%  +         ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \Psi_y  + (\partial_x v + \partial_y u) \, \partial_x \Psi_y \big ) \, \phi_i \; dx \, dy 
% $$
%
% This is a vector.
%
%
% $$
% \langle  \delta^2_{A_i\,A_j} F^x | \Psi_x \rangle + \langle  \delta^2_{A_i\,A_j} F^y | \Psi_y \rangle 
% = -\int \, \partial^2_{AA} \eta \,  \, h \, 
%   \big  (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \Psi_x  + (\partial_y u + \partial_x v) \, \partial_y \Psi_x 
%  +         ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \Psi_y  + (\partial_x v + \partial_y u) \, \partial_x \Psi_y \big ) \, \phi_i \, \phi_j \; dx \, dy 
% $$
%
% This is a matrix
%
%%


ndim=2;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psixnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psiynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

CtrlVar.EffectiveViscosity.CalculateDerivatives=false;


% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);
   

    hint=hnod*fun;
    nint=nnod*fun;
    %AGlenInt=AGlennod*fun;
    %AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    [AGlenInt,dAeffdA]=SmoothFloor(AGlennod*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);


    exx=zeros(MUA.Nele,1); exy=zeros(MUA.Nele,1);  eyy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod));


        dlxdx=dlxdx+Deriv(:,1,Inod).*Psixnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*Psiynod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*Psixnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*Psiynod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);

    [~,~,~,detadA]=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);
    detadA = detadA.*dAeffdA ;


    for Inod=1:MUA.nod
        T(:,Inod)=T(:,Inod)...
            +detadA.*hint.*((4*exx+2*eyy).*dlxdx+2*exy.*dlxdy+(4*eyy+2*exx).*dlydy+2*exy.*dlydx).*fun(Inod).*detJw;

    end
end

dIdA=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdA=dIdA+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


%% test the gradient, do the test in linear space, so do it ahead of eventual conversion to log space
if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdA);
end


%% conversion to leg space
if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    dIdA=log(10)*F.AGlen.*dIdA;
end

%% sometimes I modify the gradient to make it a L^2 or H^1 gradient instead of the l^2 gradient that I have just calculated.
% dIdA=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.A,dIdA);
% I now do this once and for all in Regularisation/Misfit







end



function   FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdA)

iNode=randi(MUA.Nnodes);

A0=F.AGlen;

deltaA=1e-3*abs(F.AGlen(iNode));

F.AGlen=A0; 
F.AGlen(iNode)=F.AGlen(iNode)-deltaA;

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=true;

Ruv_minus = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_minus = Ruv_minus.' * [Psi_x; Psi_y];

F.AGlen=A0; 
F.AGlen(iNode)=F.AGlen(iNode)+deltaA;
[Ruv_plus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_plus = Ruv_plus.' * [Psi_x; Psi_y];

F.AGlen=A0; 

dIdA_FD = (b_plus - b_minus)/(2*deltaA);   % length 2*Nnodes: top half -> column of F^{qq}_{uu}, bottom half -> column of F^{qq}_{vu}

%% there is a sign mistake which I have started to carry through the code (must correct this properly one day)
% This sign mistake has been corrected and all OK now (21/09/2026)


%[dIdA_FD dIdA(iNode)]
Diff=norm(dIdA(iNode) - dIdA_FD)/(abs(dIdA(iNode))+eps);
fprintf("dIdAq: normalized norm of difference between Direct-Adjoint and FD for node %i is %g \n",iNode,Diff)







end






