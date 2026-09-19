


function dIdC=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%% Calculates the vector quantity:
%
%
% $$ \langle  \delta_{C_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{C_i} F^y \phi_i| \Psi_y \rangle $$
%
% We expand all variables in the same FE basis 
%
% $$ \{ \phi_i \vert i=1 \ldots n \} $$
%
% and arrive at
%
%
% $$ \langle  \partial_C F^x \phi_i | \phi_j \rangle \, \Psi_{x,j} + \langle  \partial_{C} F^y \phi_i| \phi_j \rangle \, \Psi_{y,h}  $$
%
% We do not calculated the matrices 
%
% $$ \langle  \partial_C F^x \phi_i | \phi_j \rangle  $$
% 
% $$ \langle  \partial_{C} F^y \phi_i| \phi_j \rangle  $$
%
% but assembly the vector directly without the need to store these matrices.
%
%
%% Background: Formulation of the momentum equations, using the SSA approximation.
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
% $$ \rho \partial_t h + \partial_x ( \rho h u) + \partial_y (\rho h v) = \rho a $$
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
% The $uv$ forward problem is solving the equations above for $u$ and $v$. 
%
% The inverse problem is selecting one or more of $A$, $B$ and $C$ given measurements of $u$ and $v$, while respecting any
% other direct/prior information about $A$, $B$ and $C$. 
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
%%
%%
ndim=2;

% C0=CtrlVar.Czero;
% u0=CtrlVar.SpeedZero;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if CtrlVar.SlidingLaw=="Budd"
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    qnod=mnod*0 ;  % just to avoid asking this again within a loop
end

if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    muknod=mnod*0 ;
end

if ~isempty(F.V0)
    V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    V0nod=mnod*0 ;
end

Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hfnod=F.rhow*(Snod-Bnod)./rhonod;

Psi_x_nod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_y_nod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);


% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    
    
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    Cint=Cnod*fun; 
    Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    qint=qnod*fun;
    mukint=muknod*fun;
    V0int=V0nod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;
    Psi_x_int=Psi_x_nod*fun;
    Psi_y_int=Psi_y_nod*fun;
    hfint=hfnod*fun;
    %hfint=(Sint-Bint)*F.rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    %
    % dF/dC=dtaux/dC uAdjoint + dtauy/dC vAdjoint
    %
    % dtaux/dC= He u * dbeta2/dC
    %
    % beta2= (C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    %

    % setting this CtrlVar field to true ensures that BasalDrag.m returns the (point) derivative
    CtrlVar.Inverse.dFuvdClambda=true;
    Ctemp= ...
        BasalDrag(CtrlVar,MUA,Heint,[],hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,[],[],[],[],[],[],[],[],qint,F.g,mukint,V0int);
    CtrlVar.Inverse.dFuvdClambda=false;
    

    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod


        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*Psi_x_int+vint.*Psi_y_int).*fun(Inod).*detJw;  

    end
end

dIdC=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdC=dIdC+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


%% test the gradient, do the test in linear space, so do it ahead of eventual conversion to log space
if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdC);
end


%% conversion to log space
% change of variables should be done on nodal values!
% 
if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        dIdC=log(10)*F.C.*dIdC;
end

%% sometimes I modify the gradient to make it a L^2 or H^1 gradient instead of the l^2 gradient that I have just calculated.
% dIdC=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.C,dIdC);
% I now do this once and for all in Misfit/Regularisation 



end





function   FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdC)

iNode=randi(MUA.Nnodes);

C0=F.C;

deltaC=1e-3*abs(F.C(iNode));

F.C=C0; 
F.C(iNode)=F.C(iNode)-deltaC;

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=true;

Ruv_minus = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_minus = Ruv_minus.' * [Psi_x; Psi_y];

F.C=C0; 
F.C(iNode)=F.C(iNode)+deltaC;
[Ruv_plus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_plus = Ruv_plus.' * [Psi_x; Psi_y];

F.C=C0; 

dIdC_FD = (b_plus - b_minus)/(2*deltaC);   % length 2*Nnodes: top half -> column of F^{qq}_{uu}, bottom half -> column of F^{qq}_{vu}

%% there is a sign mistake which I have started to carry through the code (must correct this properly one day)
%dIdC_FD=-dIdC_FD ; % need to correct for wrong sign...


%[dIdC_FD dIdC(iNode)]
Diff=norm(dIdC(iNode) - dIdC_FD)/(abs(dIdC(iNode))+eps);
fprintf("dIdCq: normalized norm of difference between Direct-Adjoint and FD for node %i is %g \n",iNode,Diff)







end






