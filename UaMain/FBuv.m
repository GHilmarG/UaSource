

function [KFBu,KFBv]=FBuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)  %#ok<INUSD>

%% Builds the mixed Hessian term
%
% $$ \mathcal{F}^{B q} = \langle \Psi , \delta^2_{B q} \mathcal{F} \rangle , \qquad q=(u,v) $$
%
% i.e. the mixed second-order derivative of the forward model with respect to the bedrock and the velocities,
% contracted with the adjoint variables. This is the B analogue of FAuv.m and FCuv.m
%
% Two matrices are returned, each $n \times n$, with the row index corresponding to the bedrock node and the column
% index to the velocity node:
%
% $$ (K_{Bu})_{ji} = \frac{\partial^2 }{\partial B_j \, \partial u_i}
%        \Big( \langle \mathcal{F}^x , \Psi_x \rangle + \langle \mathcal{F}^y , \Psi_y \rangle \Big) $$
%
%% Why this is much simpler than FBB
%
% Write the contracted residual density as in FBB.m,
%
% $$ W = h \, V + t_{bx}\Psi_x + t_{by}\Psi_y
%      + g \cos\alpha \, (\rho h - \rho_o d) \, \Gamma
%      - \rho g \sin\alpha \, h \, \Psi_x
%      - \frac{1}{2} g \cos\alpha \, (\rho h^2 - \rho_o d^2) \, P $$
%
% Of the quantities through which B enters, only two carry any velocity dependence: the viscous term $V$ and the
% basal traction. The bed slope terms and the draft $d$ do not depend on $u$ or $v$ at all. Differentiating the
% first variation of $W$ with respect to $u_i$ therefore leaves only two terms,
%
% $$ \frac{\partial^2 W}{\partial B_j \, \partial u_i}
%      = h'_j \phi_j \, \frac{\partial V}{\partial u_i}
%      + \mu_j \phi_j \, \frac{\partial T_1}{\partial u_i} $$
%
% where $h'_j$ is the nodal derivative of the thickness from dGeometrydB.m, $\mu_j = h'_j + \rho_o/\rho$ is the
% flotation measure, and
%
% $$ T_1 = \frac{\partial t_{bx}}{\partial h}\Psi_x + \frac{\partial t_{by}}{\partial h}\Psi_y $$
%
% The $\mu_j$ arises because the drag depends on $h$ and $h_f$ only through $\Delta h = h - h_f$, so the two routes
% combine exactly as they do in the gradient. Note that no second derivative of the geometrical closure is needed
% here: $b''$ and $h''$ appear in FBB.m but not in this term.
%
%% The two velocity derivatives
%
% With $E$ the second output of EffectiveViscositySSTREAM.m,
%
% $$ \frac{\partial \eta}{\partial u_i} = E \Big[ (2\partial_x u + \partial_y v)\partial_x \phi_i
%     + \dot\epsilon_{xy} \, \partial_y \phi_i \Big] $$
%
% which is the quantity called Deu in uvMatrixAssemblySSTREAM.m, and
%
% $$ \frac{\partial V}{\partial u_i} = \frac{\partial \eta}{\partial u_i} \, \frac{V}{\eta}
%   + \eta \Big[ 4 \partial_x \phi_i \, \partial_x \Psi_x + \partial_y \phi_i \, \partial_y \Psi_x
%     + 2 \partial_x \phi_i \, \partial_y \Psi_y + \partial_y \phi_i \, \partial_x \Psi_y \Big] $$
%
% The mixed derivatives of the basal traction come from WeertmanSecondOrderDerivatives.m
%
%% Limitation
%
% The mixed thickness-velocity derivatives of the basal traction are currently available for the Weertman sliding law
% only. Everything else in this function is general.
%
%  see also: FAuv.m, FCuv.m, Hess_qp.m, FBB.m, dGeometrydB.m, WeertmanSecondOrderDerivatives.m
%
%%

narginchk(7,7)
nargoutchk(2,2)

KFBu=[] ; KFBv=[] ;

if ~contains(CtrlVar.Inverse.InvertFor,"-B-")
    % Returning empty matrices here allows Hess_qp.m to build its block matrices by concatenation, with inactive
    % fields dropping out automatically. This mirrors FAuv.m and FCuv.m
    return
end

if CtrlVar.IncludeMelangeModelPhysics
    error("FBuv:MelangeNotImplemented","FBuv is not implemented for melange model physics.")
end

if CtrlVar.Hh0~=0
    error("FBuv:NonZeroHh0",...
        "A symmetrical Heaviside function is assumed, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

ndim=2;
nNodes=MUA.Nnodes ;

%% nodal derivative of the thickness with respect to B

[~,dhdB]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow);

%% nodal values gathered onto the elements

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_xnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_ynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);

dhdBnod=reshape(dhdB(MUA.connectivity,1),MUA.Nele,MUA.nod);

dFBu=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFBv=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end

CtrlVar.EffectiveViscosity.CalculateDerivatives=true;   % E is needed below

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);

    Dx=reshape(Deriv(:,1,:),MUA.Nele,MUA.nod);
    Dy=reshape(Deriv(:,2,:),MUA.Nele,MUA.nod);

    hint=hnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;

    uint=unod*fun;
    vint=vnod*fun;

    mint=mnod*fun;
    nint=nnod*fun;

    Cint=SmoothFloor(Cnod*fun,CtrlVar.Cmin,CtrlVar.CminWidth);
    AGlenInt=SmoothFloor(AGlennod*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

    Psi_x_int=Psi_xnod*fun;
    Psi_y_int=Psi_ynod*fun;

    hfint = F.rhow*Hint./rhoint;
    Dhint = hint-hfint ;

    exx=sum(Dx.*unod,2);
    eyy=sum(Dy.*vnod,2);
    exy=0.5*(sum(Dx.*vnod,2)+sum(Dy.*unod,2));

    dlxdx=sum(Dx.*Psi_xnod,2);  dlxdy=sum(Dy.*Psi_xnod,2);
    dlydx=sum(Dx.*Psi_ynod,2);  dlydy=sum(Dy.*Psi_ynod,2);

    [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);

    % mixed thickness-velocity derivatives of the basal traction
    [~,~,d2txdhdu,d2txdhdv,d2tydhdu,d2tydhdv]=...
        WeertmanSecondOrderDerivatives(CtrlVar,Dhint,uint,vint,Cint,mint);

    detJw=detJ*MUA.weights(Iint);

    % the viscous bracket, without the viscosity
    Vbracket = (4*exx+2*eyy).*dlxdx + 2*exy.*dlxdy + (4*eyy+2*exx).*dlydy + 2*exy.*dlydx ;

    % the flotation measure, per node
    kap = F.rhow./rhoint ;

    for Jnod=1:MUA.nod       % bedrock index

        hp = dhdBnod(:,Jnod) ;
        mu = hp + kap ;

        for Inod=1:MUA.nod   % velocity index

            % d eta / d u_i   and   d eta / d v_i   (Deu and Dev of uvMatrixAssemblySSTREAM.m)
            detadu = Eint.*( (2*exx+eyy).*Dx(:,Inod) + exy.*Dy(:,Inod) ) ;
            detadv = Eint.*( (2*eyy+exx).*Dy(:,Inod) + exy.*Dx(:,Inod) ) ;

            % d V / d u_i  and  d V / d v_i
            dVdu = detadu.*Vbracket ...
                 + etaint.*( 4*Dx(:,Inod).*dlxdx + Dy(:,Inod).*dlxdy ...
                           + 2*Dx(:,Inod).*dlydy + Dy(:,Inod).*dlydx ) ;

            dVdv = detadv.*Vbracket ...
                 + etaint.*( 2*Dy(:,Inod).*dlxdx + Dx(:,Inod).*dlxdy ...
                           + 4*Dy(:,Inod).*dlydy + Dx(:,Inod).*dlydx ) ;

            % d T1 / d u_i  and  d T1 / d v_i
            dT1du = d2txdhdu.*Psi_x_int + d2tydhdu.*Psi_y_int ;
            dT1dv = d2txdhdv.*Psi_x_int + d2tydhdv.*Psi_y_int ;

            dFBu(:,Jnod,Inod)=dFBu(:,Jnod,Inod) ...
                + ( hp.*dVdu + mu.*dT1du.*fun(Inod) ).*fun(Jnod).*detJw ;

            dFBv(:,Jnod,Inod)=dFBv(:,Jnod,Inod) ...
                + ( hp.*dVdv + mu.*dT1dv.*fun(Inod) ).*fun(Jnod).*detJw ;

        end
    end
end

%% assemble

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Uval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Vval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Jnod=1:MUA.nod
    for Inod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);   % rows : bedrock nodes
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);   % cols : velocity nodes
        Uval(istak+1:istak+MUA.Nele)=dFBu(:,Jnod,Inod);
        Vval(istak+1:istak+MUA.Nele)=dFBv(:,Jnod,Inod);
        istak=istak+MUA.Nele;
    end
end

KFBu=sparse(Iind,Jind,Uval,nNodes,nNodes);
KFBv=sparse(Iind,Jind,Vval,nNodes,nNodes);

end
