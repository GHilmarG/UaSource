

function KFBB=FBB(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)  %#ok<INUSD>

%% Builds the Hessian term
%
% $$ \mathcal{F}^{BB} = \langle \Psi , \delta^2_{BB} \mathcal{F} \rangle $$
%
% i.e. the second-order derivative of the forward model with respect to the bedrock, contracted with the adjoint
% variables. This is the B analogue of FAA.m and FCC.m
%
%% Structure of the calculation
%
% Write the contracted residual density as
%
% $$ W = h \, V + t_{bx}\Psi_x + t_{by}\Psi_y
%      + g \cos\alpha \, (\rho h - \rho_o d) \, \Gamma
%      - \rho g \sin\alpha \, h \, \Psi_x
%      - \frac{1}{2} g \cos\alpha \, (\rho h^2 - \rho_o d^2) \, P $$
%
% with
%
% $$ V = \eta \Big[ (4 \partial_x u + 2\partial_y v)\partial_x \Psi_x + (\partial_y u + \partial_x v)\partial_y \Psi_x
%        + (4 \partial_y v + 2\partial_x u)\partial_y \Psi_y + (\partial_y u + \partial_x v)\partial_x \Psi_y \Big] $$
%
% $$ \Gamma = \partial_x b \, \Psi_x + \partial_y b \, \Psi_y , \qquad P = \partial_x \Psi_x + \partial_y \Psi_y $$
%
% Neither $V$ nor $P$ depends on B.
%
%% Only four quantities carry the B dependence
%
% At an integration point, B enters $W$ only through
%
% $$ h , \qquad H=S-B , \qquad \partial_x b , \qquad \partial_y b $$
%
% since the flotation measure is $\Delta h = h - \kappa H$ with $\kappa=\rho_o/\rho$, and the draft is $d=d(h,H)$.
% The derivatives of these four with respect to the nodal values of B are
%
% $$ \frac{\partial h}{\partial B_j} = h'_j \phi_j , \qquad
%    \frac{\partial H}{\partial B_j} = -\phi_j , \qquad
%    \frac{\partial (\partial_x b)}{\partial B_j} = b'_j \, \partial_x \phi_j $$
%
% and the only non-zero second derivatives are the diagonal ones
%
% $$ \frac{\partial^2 h}{\partial B_j \partial B_l} = h''_j \, \delta_{jl} \, \phi_j , \qquad
%    \frac{\partial^2 (\partial_x b)}{\partial B_j \partial B_l} = b''_j \, \delta_{jl} \, \partial_x \phi_j , \qquad
%    \frac{\partial^2 H}{\partial B_j \partial B_l} = 0 $$
%
% with $h'$, $b'$, $h''$ and $b''$ supplied by dGeometrydB.m. They are nodal scalars, because the geometrical
% closure is applied independently at each node, which is why the second-derivative contributions are diagonal.
%
% Collecting, and writing $X^{\alpha}_j$ for the four first derivatives above,
%
% $$ \frac{\partial^2 W}{\partial B_j \partial B_l} = \sum_{\alpha\beta} W_{\alpha\beta} X^\alpha_j X^\beta_l
%    \;+\; \delta_{jl} \left ( W_h \, h''_j \phi_j + W_{b_x} b''_j \partial_x \phi_j + W_{b_y} b''_j \partial_y \phi_j \right ) $$
%
% Since $W$ is LINEAR in $\partial_x b$ and $\partial_y b$, the blocks $W_{b_x b_x}$, $W_{b_x b_y}$ and
% $W_{b_y b_y}$ all vanish, which removes six of the ten possible second-derivative combinations.
%
%% Derivatives of the draft
%
% With $\Theta = \tilde{\delta}(\Delta h) ( H^{+} - \rho h/\rho_o )$, $\Theta' = \tilde{\delta}'(\Delta h) ( H^{+} - \rho h/\rho_o )$,
% and $H^{+}=\mathcal{H}(H) H$ so that $(H^{+})' = \tilde{\delta}(H) H + \mathcal{H}(H)$ and
% $(H^{+})'' = \tilde{\delta}'(H) H + 2\tilde{\delta}(H)$,
%
% $$ d_h = \Theta + (1-\mathcal{G}) \frac{\rho}{\rho_o} , \qquad  d_H = -\kappa \Theta + \mathcal{G} (H^{+})' $$
%
% $$ d_{hh} = \Theta' - 2 \tilde{\delta}(\Delta h) \frac{\rho}{\rho_o} , \qquad
%    d_{hH} = -\kappa \Theta' + \tilde{\delta}(\Delta h) (H^{+})' + \kappa \tilde{\delta}(\Delta h) \frac{\rho}{\rho_o} $$
%
% $$ d_{HH} = \kappa^2 \Theta' - 2 \kappa \tilde{\delta}(\Delta h) (H^{+})' + \mathcal{G} (H^{+})'' $$
%
% Here $\tilde{\delta}$ is the derivative of the SMOOTHED Heaviside function, not a Dirac delta, so $\tilde{\delta}'$
% is an ordinary classical derivative, obtained from $\tilde{\delta}'=2k\tilde{\delta}(1-2\mathcal{H})$.
%
%% Limitation
%
% The second derivative of the basal traction with respect to h is obtained from
% WeertmanSecondOrderDerivatives.m, which currently supports the Weertman sliding law only. Everything else in this
% function is general.
%
%  see also: FAA.m, FCC.m, Fpp.m, dGeometrydB.m, WeertmanSecondOrderDerivatives.m, dIdBqGeneral.m
%
%%

narginchk(7,7)
nargoutchk(1,1)

if CtrlVar.IncludeMelangeModelPhysics
    error("FBB:MelangeNotImplemented","FBB is not implemented for melange model physics.")
end

if CtrlVar.Hh0~=0
    error("FBB:NonZeroHh0",...
        "A symmetrical Heaviside function is assumed, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

ndim=2;
ca=cos(F.alpha) ; sa=sin(F.alpha) ;
nNodes=MUA.Nnodes ;
kH=CtrlVar.kH ;

%% nodal derivatives of the geometry with respect to B

[dbdB,dhdB,~,~,d2bdB2,d2hdB2]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow);

%% nodal values gathered onto the elements

snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
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

dbdBnod  =reshape(dbdB(MUA.connectivity,1),MUA.Nele,MUA.nod);
dhdBnod  =reshape(dhdB(MUA.connectivity,1),MUA.Nele,MUA.nod);
d2bdB2nod=reshape(d2bdB2(MUA.connectivity,1),MUA.Nele,MUA.nod);
d2hdB2nod=reshape(d2hdB2(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F.q)   ; qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod)     ; else ; qnod=[]   ; end
if ~isempty(F.muk) ; muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod) ; else ; muknod=[] ; end
if ~isempty(F.V0)  ; V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod)   ; else ; V0nod=[]  ; end

d2FBB=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end

CtrlVar.BasalDrag.CalculateDerivatives=true;

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

    % C and AGlen are coefficients here, B is the differentiation variable
    Cint=SmoothFloor(Cnod*fun,CtrlVar.Cmin,CtrlVar.CminWidth);
    AGlenInt=SmoothFloor(AGlennod*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

    if ~isempty(qnod)   ; qint=qnod*fun     ; else ; qint=[]   ; end
    if ~isempty(muknod) ; mukint=muknod*fun ; else ; mukint=[] ; end
    if ~isempty(V0nod)  ; V0int=V0nod*fun   ; else ; V0int=[]  ; end

    Psi_x_int=Psi_xnod*fun;
    Psi_y_int=Psi_ynod*fun;

    % masks
    hfint=F.rhow*Hint./rhoint;
    Dhint    = hint-hfint ;
    Heint    = HeavisideApprox(kH,Dhint,CtrlVar.Hh0);
    deltaint = DiracDelta(kH,Dhint,CtrlVar.Hh0);
    dDelta   = 2*kH*deltaint.*(1-2*Heint) ;                 % derivative of the regularised delta

    HeHint    = HeavisideApprox(kH,Hint,CtrlVar.Hh0);
    deltaHint = DiracDelta(kH,Hint,CtrlVar.Hh0);
    dDeltaH   = 2*kH*deltaHint.*(1-2*HeHint) ;
    Hposint   = HeHint.*Hint;

    dint = (1-Heint).*rhoint.*hint/F.rhow + Heint.*Hposint ;

    kap = F.rhow./rhoint ;      % rho_o/rho
    rr  = rhoint/F.rhow ;       % rho/rho_o

    % gradients
    dsdx=sum(Dx.*snod,2);   dsdy=sum(Dy.*snod,2);
    dhdx=sum(Dx.*hnod,2);   dhdy=sum(Dy.*hnod,2);
    dbdx=dsdx-dhdx;         dbdy=dsdy-dhdy;

    exx=sum(Dx.*unod,2);
    eyy=sum(Dy.*vnod,2);
    exy=0.5*(sum(Dx.*vnod,2)+sum(Dy.*unod,2));

    dlxdx=sum(Dx.*Psi_xnod,2);  dlxdy=sum(Dy.*Psi_xnod,2);
    dlydx=sum(Dx.*Psi_ynod,2);  dlydy=sum(Dy.*Psi_ynod,2);

    etaint=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);

    [~,~,~,~,~,~,dtaubxdh,dtaubydh] = ...
        BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,...
        [],[],[],[],[],[],[],[],qint,F.g,mukint,V0int);

    [d2taubxdh2,d2taubydh2]=WeertmanSecondOrderDerivatives(CtrlVar,Dhint,uint,vint,Cint,mint);

    detJw=detJ*MUA.weights(Iint);

    %% derivatives of the draft with respect to h and H

    HmRh  = Hposint - rhoint.*hint/F.rhow ;     % H^+ - rho h/rho_o
    Theta = deltaint.*HmRh ;
    ThetaP= dDelta.*HmRh ;

    Hp1 = deltaHint.*Hint + HeHint ;            % (H^+)'
    Hp2 = dDeltaH.*Hint + 2*deltaHint ;         % (H^+)''

    d_h  = Theta + (1-Heint).*rr ;
    d_H  = -kap.*Theta + Heint.*Hp1 ;

    d_hh = ThetaP - 2*deltaint.*rr ;
    d_hH = -kap.*ThetaP + deltaint.*Hp1 + kap.*deltaint.*rr ;
    d_HH = kap.^2.*ThetaP - 2*kap.*deltaint.*Hp1 + Heint.*Hp2 ;

    %% contractions with the adjoint fields

    V  = etaint.*( (4*exx+2*eyy).*dlxdx + 2*exy.*dlxdy + (4*eyy+2*exx).*dlydy + 2*exy.*dlydx ) ;
    Gm = dbdx.*Psi_x_int + dbdy.*Psi_y_int ;
    Pr = dlxdx + dlydy ;
    T1 = dtaubxdh.*Psi_x_int  + dtaubydh.*Psi_y_int ;
    T2 = d2taubxdh2.*Psi_x_int + d2taubydh2.*Psi_y_int ;

    gc = F.g*ca ;  gs = F.g*sa ;
    Buoy = rhoint.*hint - F.rhow*dint ;

    %% partial derivatives of W

    W_h  = V + T1 + gc*(rhoint - F.rhow*d_h).*Gm - gs*rhoint.*Psi_x_int ...
           - gc*(rhoint.*hint - F.rhow*dint.*d_h).*Pr ;
    % Note: W_H, the first derivative of W with respect to H, is deliberately absent. It would only ever multiply
    % the second derivative of H with respect to B, and H=S-B is linear in B, so that second derivative is exactly
    % zero. Contrast W_h and W_bx/W_by, which do appear below, because h and the bed slope depend on B through the
    % non-linear geometrical closure and therefore have non-zero second derivatives.
    % W_H  = -kap.*T1 - gc*F.rhow*d_H.*Gm + gc*F.rhow*dint.*d_H.*Pr ;

    W_bx = gc*Buoy.*Psi_x_int ;
    W_by = gc*Buoy.*Psi_y_int ;

    W_hh = T2 - gc*F.rhow*d_hh.*Gm - gc*( rhoint - F.rhow*(d_h.^2 + dint.*d_hh) ).*Pr ;
    W_hH = -kap.*T2 - gc*F.rhow*d_hH.*Gm + gc*F.rhow*( d_h.*d_H + dint.*d_hH ).*Pr ;
    W_HH = kap.^2.*T2 - gc*F.rhow*d_HH.*Gm + gc*F.rhow*( d_H.^2 + dint.*d_HH ).*Pr ;

    W_hbx = gc*(rhoint - F.rhow*d_h).*Psi_x_int ;
    W_hby = gc*(rhoint - F.rhow*d_h).*Psi_y_int ;
    W_Hbx = -gc*F.rhow*d_H.*Psi_x_int ;
    W_Hby = -gc*F.rhow*d_H.*Psi_y_int ;

    %% assemble

    for Jnod=1:MUA.nod

        Xh_J  = dhdBnod(:,Jnod).*fun(Jnod) ;
        XH_J  = -fun(Jnod) ;
        Xbx_J = dbdBnod(:,Jnod).*Dx(:,Jnod) ;
        Xby_J = dbdBnod(:,Jnod).*Dy(:,Jnod) ;

        for Lnod=1:MUA.nod

            Xh_L  = dhdBnod(:,Lnod).*fun(Lnod) ;
            XH_L  = -fun(Lnod) ;
            Xbx_L = dbdBnod(:,Lnod).*Dx(:,Lnod) ;
            Xby_L = dbdBnod(:,Lnod).*Dy(:,Lnod) ;

            T =   W_hh.*Xh_J.*Xh_L ...
                + W_hH.*( Xh_J.*XH_L + XH_J.*Xh_L ) ...
                + W_HH.*( XH_J.*XH_L ) ...
                + W_hbx.*( Xh_J.*Xbx_L + Xbx_J.*Xh_L ) ...
                + W_hby.*( Xh_J.*Xby_L + Xby_J.*Xh_L ) ...
                + W_Hbx.*( XH_J.*Xbx_L + Xbx_J.*XH_L ) ...
                + W_Hby.*( XH_J.*Xby_L + Xby_J.*XH_L ) ;

            if Jnod==Lnod
                % the only second derivatives of h and of the bed slope are diagonal, because the geometrical
                % closure is applied independently at each node
                T = T + W_h.*d2hdB2nod(:,Jnod).*fun(Jnod) ...
                      + W_bx.*d2bdB2nod(:,Jnod).*Dx(:,Jnod) ...
                      + W_by.*d2bdB2nod(:,Jnod).*Dy(:,Jnod) ;
            end

            d2FBB(:,Jnod,Lnod)=d2FBB(:,Jnod,Lnod) + T.*detJw ;

        end
    end
end

%% assemble into a sparse n x n matrix

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Jnod=1:MUA.nod
    for Lnod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Lnod);
        Xval(istak+1:istak+MUA.Nele)=d2FBB(:,Jnod,Lnod);
        istak=istak+MUA.Nele;
    end
end

KFBB=sparse(Iind,Jind,Xval,nNodes,nNodes);

KFBB=(KFBB+KFBB.')/2 ;   % symmetric by construction, this removes round-off asymmetry

end
