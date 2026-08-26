


function dIdC=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%% Calculates the vector quantity:
%
%
% $$ \langle  \delta_{C_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{C_i} F^y \phi_i| \Psi_y \rangle $$
%
%
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
%
%%
ndim=2;

C0=CtrlVar.Czero;
u0=CtrlVar.SpeedZero;

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

    % this is a temporary change which only works for Weertman sliding law, 
     % speed=sqrt(uint.*uint+vint.*vint+u0^2);
     % Um=speed.^(1./mint-1) ;
     % Ctemp =  Heint.*    (1./mint).*(Cint+C0).^(-1./mint-1)   .*Um;  % Um=speed.^(1./m-1) ; This is the same Ctemp as returned by BasalDrag for Weertman sliding law, 
    % but the BasalDrag function is more general and will return this quantity for various other sliding laws as well 
    %%


    % Note: I include the u and v in the adjoint calculation itself below, so I just need the
    % derivative without the u and the v. Therefore 
    %


    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod


        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*Psi_x_int+vint.*Psi_y_int).*fun(Inod).*detJw;  

    end
end

dIdC=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdC=dIdC+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


% change of variables should be done on nodal values!
% I learned this the hard way by doing extensive tests on dJ/dgamma
if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        dIdC=log(10)*F.C.*dIdC;
end


% dIdC=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.C,dIdC);



end





