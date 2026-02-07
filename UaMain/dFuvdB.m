function K=dFuvdB(CtrlVar,MUA,F)

%%
%
% assembles the matrix K which is the FE form of
%
% $$d_{\mathbf{B}} \mathbf{F}_{\mathbf{v}}$$
%
% where
%
% $$\mathbf{F}_{\mathbf{v}}$$
%
% is the uv-forward model.
%
%
% $$\left[\begin{array}{cccc}
% \partial F^x_1 /\partial B_1  & \partial F^x_1 /\partial B_2  & \ldots & \partial F^x_1 /\partial B_n  \\
% \partial F^x_2 /\partial B_1  & \partial F^x_2 /\partial B_2  & \ldots & \partial F^x_2 /\partial B_n  \\
%              \vdots              &              \vdots              &  \vdots  &    \vdots              \\
% \partial F^y_1 /\partial B_1  & \partial F^y_1 /\partial B_2 & \ldots & \partial F^y_1 /\partial B_n  \\
% \partial F^y_2 /\partial B_1  & \partial F^y_2 /\partial B_2 & \ldots & \partial F^y_2 /\partial B_n  \\
%              \vdots              &              \vdots              &  \vdots  &    \vdots              \\
% \end{array}\right] $$
%
%
%
% $$ F_x=\partial_x ( h \eta ( 4 \partial_x u + 2 \partial_y v)) + \partial_y ( h \eta (\partial_y u + \partial_x v) ) - t_x
% -   \frac{1}{2} g \cos(\alpha) \partial_x (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+})\, \partial_x B + \rho g \sin(\alpha) \, h =0
% $$
%
% $$ F_y = \partial_y ( h \eta ( 4 \partial_y v + 2 \partial_x u)) + \partial_x ( h \eta (\partial_x v + \partial_y y) ) - t_y
% -   \frac{1}{2} g \cos(\alpha) \partial_y (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+})\, \partial_y B =0
% $$
%
% or
%
% $$ F_x=\partial_x ( h \eta ( 4 \partial_x u + 2 \partial_y v)) + \partial_y ( h \eta (\partial_y u + \partial_x v) ) - t_x
% -   \frac{1}{2} g \cos(\alpha) \partial_x (\rho h^2 -  \rho_o d^2)- g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_x B+ \rho g \sin(\alpha) \, h  =0
% $$
%
% $$ F_y = \partial_y ( h \eta ( 4 \partial_y v + 2 \partial_x u)) + \partial_x ( h \eta (\partial_x v + \partial_y y) ) - t_y
% - \frac{1}{2} g \cos(\alpha) \partial_y (\rho h^2 -  \rho_o d^2)- g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =0
% $$
%
%
% with
%
% $$ t_x = t_x(C,u,v) $$
%
% $$ t_y = t_y(C,u,v) $$
%
% and
%
% $$ \eta=\eta(A,u,v) $$
%
%
%
% $$ \partial F^x_i/\partial A_j =  \langle h (\partial \eta/\partial A)  ( 4 \partial_x u + 2 \partial_y v))\, \phi_j | \partial_x \phi_i \rangle
% + \langle h (\partial \eta/\partial A) (\partial_y u + \partial_x v)\, \phi_j | \partial_y \phi_i \rangle $$
%
% $$ \partial F^y_i/\partial A_j =  \langle h (\partial \eta/\partial A)  ( 4 \partial_y v + 2 \partial_x u))\, \phi_j | \partial_y \phi_i \rangle
% + \langle h (\partial \eta/\partial A) (\partial_x v + \partial_y u)\, \phi_j | \partial_x \phi_i \rangle $$
%
% Currently only done for grounded area where $B=b$, but once this works for grounded area, will do for the more general
% case...
%
% $$b=\mathcal{G} B + (1-\mathcal{G}) (S-\rho h/\rho_0) $$
%
% The above equation for $b$ contains $h$ on the right hand side as well through $h=s-b$. This is OK if we consider $h$ as
% given. Here we consider the upper surface $s$ as given (i.e. as a input) and therefore solve for $b$, arriving at
%
% $$ b=\frac{1}{1-(1-\mathcal{G}) \rho/\rho_o} \, ( \mathcal{G} B + (1-\mathcal{G}) ( S - \rho s/\rho_o ) ) $$
%
% However, this expression still contains a dependency on $b$ on the right-hand side because $\mathcal{G}$ depends on $h$ and
% therefore on $b$. This is, thus, an implicit function for $b$. To calculate the (total) derivative $db/dB$ we can use implicit differentiation.
%
%
% $$h=s-b$$
%
% $$H=S-B$$
%
% $$h_f=\rho_o H/\rho$$
%
% $$\mathcal{G}=\mathcal{H}(h-h_f)$$
%
% $$H^{+}=\mathcal{H}(H) \, H $$
%
% $$
% d=\mathcal{H}(H) \, (S-b)
% = \mathcal{H}(h_f-h) \rho h /\rho_o + \mathcal{H}(h-h_f) \; \mathcal{H}(H) \, H
% = \mathcal{H}(h_f-h) \rho h /\rho_o + \mathcal{H}(h-h_f) \; H^{+}
% = (1-\mathcal{G}) \rho h /\rho_o + \mathcal{G} \; H^{+}
% $$
%
% $$ \tau = \mathcal{G} \, \beta^2 v  $$
%
%
% Derivatives for grounded ice only, i.e. where $\mathcal{G}=1$ and $\mathrm{d}\mathcal{G}/\mathrm{d}B=0$ with $s$ and $S$ given:
%
% $$\mathrm{d}s/\mathrm{d}B=\mathrm{d}S/\mathrm{d}B=0$$
%
% $$\mathrm{d}b/\mathrm{d}B=1$$
%
% $$\mathrm{d}h/\mathrm{d}B=\mathrm{d}s/\mathrm{B}-\mathrm{d}b/\mathrm{d}B = 0-1=-1$$
%
%
% $$\mathrm{d}H/\mathrm{d}B=-1 $$
%
% $$\mathrm{d} h_f/\mathrm{d}B=(\rho_o/\rho) \; \mathrm{d}H/\mathrm{d}B = -\frac{\rho_0}{\rho} $$
%
% $$ \mathrm{d}\mathcal{G}/\mathrm{d}B = 0 $$
%
% $$ \mathrm{d} H^{+} / \mathrm{d}B = \delta(H) \; \mathrm{d}H/\mathrm{d}B \; H + \mathcal{H}(H) \; \mathrm{d}H/\mathrm{d}B = -\delta(H)\, H - \mathcal{H}(H) $$
%
% $$ \mathrm{d} d/\mathrm{d} B = \mathrm{d} H^{+}/\mathrm{d} B = -\delta(H)\, H - \mathcal{H}(H)$$
%
% $$\mathrm{d} \tau_x /\mathrm{d} B = \frac{\partial \tau_x}{\partial h } \, \frac{\mathrm{d} h}{\mathrm{d} B} = 0 \; (-1) = 0  $$
%
% FE-formulation, note how I use
%
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
%
% $$ F^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle t_x | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  + \langle \rho g \sin(\alpha) \, h  | \phi \rangle   =0
% $$
%
% $$ F^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle t_y | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
%
%
% Grounded:
%
% $$ \partial F^x_i/\partial B_j
% = \langle \eta \, (\partial h/\partial B) \, (4 \partial_x u+2 \partial_y v) \,  \phi_j | \partial_x \phi_i \rangle
% + \langle \eta \, (\partial h/\partial B) \, (\partial_y u+\partial_x v) \, \phi_j | \partial_y \phi_i \rangle
% - \langle g \cos(\alpha) \, ( \rho h \, (\partial h/\partial B) - \rho_o d \, (\partial d/\partial B) ) \,  \phi_j | \partial_x \phi_i \rangle
% + \langle g \cos(\alpha) \, \mathcal{G} \, ( \rho \,(\partial h /\partial B) - \rho_o \, (\partial H^{+}/\partial B)) \, \partial_x B \; \phi_j | \phi_i \rangle
% + \langle  g \cos(\alpha) \, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_x \phi_j | \phi_i  \rangle + \langle \rho g \sin(\alpha) \, (\partial h/\partial B) \phi_j  | \phi_i \rangle
% $$
%
% $$ \partial F^y_i/\partial B_j
% = \langle \eta \, (\partial h/\partial B) \, (4 \partial_y v+2 \partial_y u) \,  \phi_j | \partial_y \phi_i \rangle
% + \langle \eta \, (\partial h/\partial B) \, (\partial_x v+\partial_y u) \, \phi_j | \partial_x \phi_i \rangle
% - \langle g \cos(\alpha) \, ( \rho h \, (\partial h/\partial B) - \rho_o d \, (\partial d/\partial B) ) \,  \phi_j | \partial_y \phi_i \rangle
% + \langle g \cos(\alpha) \, \mathcal{G} \, ( \rho \,(\partial h /\partial B) - \rho_o \, (\partial H^{+}/\partial B)) \, \partial_y B \; \phi_j | \phi_i \rangle
% + \langle  g \cos(\alpha) \, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_y \phi_j | \phi_i  \rangle
% $$
%
%
%
% see also: duvdAFunc.m, duvdBFunc.m, duvdCFunc.m, dFuvdA.m, dFuvdB.m, dFuvdC.m, TestSensitivityMatrixCalculations.m
%
%%

%%


ndim=2;
nNodes=MUA.Nnodes ;



Anod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

ca=cos(F.alpha); sa=sin(F.alpha);

g=F.g;

dFxdB=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFydB=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end


for Iint=1:MUA.nip   % integration points

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);


    sint=snod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;

    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"

        bint=Bint ;     % ~OK, expect when grounded
        hint=sint-bint; % OK

    else    % CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;

        hint=hnod*fun;  %  I could put calculating bs from hBS in here
        bint=sint-hint; %

    end



    rhoint=rhonod*fun;
    rhoo=F.rhow;



    Aint=Anod*fun;
    Aint(Aint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    nint=nnod*fun;

    dsdx=zeros(MUA.Nele,1); dhdx=zeros(MUA.Nele,1);
    dsdy=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
    dBdx=zeros(MUA.Nele,1); dBdy=zeros(MUA.Nele,1);
    exx=zeros(MUA.Nele,1);
    eyy=zeros(MUA.Nele,1);
    exy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);

        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);

        dBdx=dBdx+Deriv(:,1,Inod).*Bnod(:,Inod);
        dBdy=dBdy+Deriv(:,2,Inod).*Bnod(:,Inod);

        exx=exx+Deriv(:,1,Inod).*ubnod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vbnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vbnod(:,Inod) + Deriv(:,2,Inod).*ubnod(:,Inod));

    end


    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"
        
        dbdx=dBdx;  % only OK if grounded
        dbdy=dBdy;
    else
        dbdx=dsdx-dhdx;   %  I could put calculating bs from hBS in here
        dbdy=dsdy-dhdy;
    end



    hfint=rhoo*Hint./rhoint;                                         % this is linear, so fine to evaluate at int in this manner
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);      % important to calculate Heint and deltaint in a consistent manner
    HEint = HeavisideApprox(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);

    dHeintdh=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    dHEintdh=-DiracDelta(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);

    G=Heint;


    Eta=EffectiveViscositySSTREAM(CtrlVar,Aint,nint,exx,eyy,exy);

    HeH=HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0);
    dHeHdH=DiracDelta(CtrlVar.kH,Hint,CtrlVar.Hh0);

    dHdB=-1;
    Hposint=HeH.*Hint;
    dHposintdB=dHeHdH.*dHdB.*Hint+HeH.*dHdB;

    dbdB=1;
    dhdb=-1 ;
    dhdB=dhdb.*dbdB;

    dhfdB=rhoo.*dHdB./rhoint;
    dGdB=0;

    dint=HEint.*rhoint.*hint/rhoo + Heint.*Hposint ;  % definition of d, applied directly at integration points, also: d=(1-G).*rhoint.Uhint./rhoo+G.*Hposint;
    dddB=dHEintdh.*dhdB.*rhoint.*hint/rhoo+HEint.*rhoint.*dhdB/rhoo+dHeintdh.*dhdB.*Hposint+Heint.*dHposintdB ;

    % taub?
    %
    % [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh,taubxo,taubyo,taubxa,taubya] = ...
    %     BasalDrag(CtrlVar,MUA,He,delta,h,B,H,rho,rhow,ub,vb,C,m,uo,vo,Co,mo,ua,va,Ca,ma,q,g,muk,V0)


    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod


            dFxdB(:,Inod,Jnod)=dFxdB(:,Inod,Jnod)+...
                (...
                + Eta.*dhdB.*(4*exx+2*eyy).*fun(Jnod).*Deriv(:,1,Inod) ...
                + Eta.*dhdB.*(2*exy).*fun(Jnod).*Deriv(:,2,Inod)...
                -g.*ca.*(rhoint.*hint.*dhdB-rhoo.*dint.*dddB).*fun(Jnod).*Deriv(:,1,Inod)...
                +g.*ca.*G.*(rhoint.*dhdB-rhoo.*dddB).*dbdx.*fun(Jnod).*fun(Inod)...
                +g.*ca.*G.*(rhoint.*hint-rhoo.*dint).*Deriv(:,1,Jnod).*fun(Inod)...
                +g.*sa.*rhoint.*dhdB.*fun(Jnod).*fun(Inod)...
                ).*detJw;

            dFydB(:,Inod,Jnod)=dFydB(:,Inod,Jnod)+...
                (...
                + Eta.*dhdB.*(4*eyy+2*exx).*fun(Jnod).*Deriv(:,2,Inod) ...
                + Eta.*dhdB.*(2*exy).*fun(Jnod).*Deriv(:,1,Inod)...
                -g.*ca.*(rhoint.*hint.*dhdB-rhoo.*dint.*dddB).*fun(Jnod).*Deriv(:,2,Inod)...
                +g.*ca.*G.*(rhoint.*dhdB-rhoo.*dddB).*dbdy.*fun(Jnod).*fun(Inod)...
                +g.*ca.*G.*(rhoint.*hint-rhoo.*dint).*Deriv(:,2,Jnod).*fun(Inod)...
                ).*detJw;




        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFxdB(:,Inod,Jnod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+nNodes; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFydB(:,Inod,Jnod);
        istak=istak+MUA.Nele;

    end
end

K=sparseUA(Iind,Jind,Xval,2*nNodes,nNodes);
end





