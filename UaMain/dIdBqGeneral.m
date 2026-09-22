


function dIdB=dIdBqGeneral(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

%% Calculates the vector quantity:
%
% $$ \langle  \delta_{B_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{B_i} F^y \phi_i| \Psi_y \rangle $$
%
% This is the general version of
%
%   dIdbq.m
%
% Unlike dIdbq.m, no assumption is made about the ice being grounded. The floating mask is itself a function of B and
% its variation with B is accounted for.
%
%% Assumptions
%
% # The upper ice surface, s, the ocean surface, S, and the densities are held fixed, and b and h are obtained from s,
%   S and B using the closure solved by Calc_bh_From_sBS.m
% # The Heaviside function is symmetrical, i.e. CtrlVar.Hh0=0 (see below).
% # All sliding laws are supported. Melange model physics is not.
%
% The first of these requires that the forward model is run with
%
%   CtrlVar.Calculate.Geometry="bs-FROM-hBS"
%
% so that the assembly uses the nodal h and b fields, rather than setting b=B at the integration points, which is only
% correct where grounded. F.b and F.h must be the converged output of Calc_bh_From_sBS.m for the current F.B. This is
% done within p2F.m at every cost-function evaluation.
%
%% The symmetrical Heaviside function
%
% With CtrlVar.Hh0=0 the smoothed Heaviside is centred on the origin,
%
% $$ \mathcal{H}(\xi)=\frac{1}{2}\left ( 1 + \tanh (k \xi) \right ) , \qquad \delta(\xi) = \mathcal{H}'(\xi) = \frac{k}{2} \, \mathrm{sech}^2 (k \xi) $$
%
% and therefore has the two properties
%
% $$ \mathcal{H}(-\xi) = 1-\mathcal{H}(\xi) , \qquad \delta(-\xi)=\delta(\xi) $$
%
% Both are used below. They mean that the two masks appearing in the draft are not independent, and that only one
% $\delta$ is ever needed. (Numerically, the second of these holds to the last bit, and the first to one ulp.)
%
%% The finite-element form of the forward model
%
% $$
% F^x_i= \left \langle  h \eta \, ( 4 \partial_x u + 2 \partial_y v) \vert  \, \partial_x \phi_i \right \rangle
% + \langle   h \eta \, (\partial_y u + \partial_x v)  \vert  \partial_y \phi_i \rangle
% + \langle t_{bx} , \phi_i \rangle
%  - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \Big\vert \partial_x \phi_i \right \rangle
% + \langle g \cos(\alpha) \, (\rho h -\rho_o d) \, \partial_x b \vert  \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% and correspondingly for $F^y_i$, without the $\sin \alpha$ term.
%
%% Dependency on B
%
% B enters in four distinct ways:
%
% # through the thickness, h, and the base, b, via the geometrical closure. These are NODAL dependencies, i.e.
%   $\partial h(x)/\partial B_j = (\partial h/\partial B)_j \, \phi_j(x)$, with the nodal derivatives supplied by
%   dGeometrydB.m
% # through the flotation thickness evaluated at the integration points, $\partial h_f/\partial B_j = -(\rho_o/\rho) \phi_j$
% # through $H=S-B$ evaluated at the integration points, $\partial H/\partial B_j = -\phi_j$
% # through $\partial_x b$, and since b is a finite-element field,
%   $\partial (\partial_x b)/\partial B_j = (\partial b/\partial B)_j \, \partial_x \phi_j$
%
%% The flotation measure
%
% Every quantity that depends on whether the ice is grounded does so through the single scalar
%
% $$ \Delta h = h - h_f $$
%
% and therefore through
%
% $$ \mu_j := \frac{\partial \Delta h}{\partial B_j} = \left ( \frac{\partial h}{\partial B} \right )_j + \frac{\rho_o}{\rho} $$
%
% where the second term comes from $\partial h_f/\partial B = -\rho_o/\rho$. This factor appears in the grounding mask,
% in the effective pressure of the sliding law, and in the draft. It is formed once per node as mu below.
%
%% Basal drag
%
% No per-sliding-law treatment is needed. Both the floating mask and the effective pressure depend on h and h_f only
% through $\Delta h$, and therefore $\partial t_{bx}/\partial h_f = -\partial t_{bx}/\partial h$ for every law. Hence
%
% $$ \frac{\partial t_{bx}}{\partial B_j} = \mu_j \, \frac{\partial t_{bx}}{\partial h} \, \phi_j $$
%
% and the quantity dtaubxdh returned by BasalDrag.m is all that is required.
%
%% The draft
%
% Using $\mathcal{H}(h_f-h)=1-\mathcal{G}$, the draft can be written with a single mask,
%
% $$ d = (1-\mathcal{G}) \, \frac{\rho h}{\rho_o} + \mathcal{G} \, H^{+} , \qquad H^{+}=\mathcal{H}(H) \, H $$
%
% B reaches d through $\Delta h$ in the mask, through the h appearing explicitly in $\rho h/\rho_o$, and through
% $H=S-B$. With
%
% $$ \Theta := \delta(\Delta h) \left ( H^{+} - \frac{\rho h}{\rho_o} \right ) $$
%
% the three routes combine to give
%
% $$ d'_j = \Theta \, \mu_j + (1-\mathcal{G}) \, \frac{\rho}{\rho_o} \, \left ( \frac{\partial h}{\partial B} \right )_j - \mathcal{G} \left ( \delta(H) \, H + \mathcal{H}(H) \right ) $$
%
% the three terms being grounding-line migration, the thinning of floating ice, and the change in the depth of the bed
% below sea level. Note that $\Theta$ is $\delta$ multiplied by the DIFFERENCE between the two candidate drafts, and so
% vanishes at flotation where the two agree. This is what keeps the draft smooth across the grounding line.
%
%  see also: dGeometrydB.m, dIdbq.m, dIdCq.m, dIdAq.m
%
%%

narginchk(7,7)

if CtrlVar.IncludeMelangeModelPhysics
    error("dIdBqGeneral:MelangeNotImplemented","Inversion for B with melange model physics is not implemented.")
end

if CtrlVar.Hh0~=0
    error("dIdBqGeneral:NonZeroHh0",...
        "This implementation assumes a symmetrical Heaviside function, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

ndim=2;
ca=cos(F.alpha) ; sa=sin(F.alpha) ;

%% nodal derivatives of the geometry with respect to B

[dbdB,dhdB]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow);

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

dbdBnod=reshape(dbdB(MUA.connectivity,1),MUA.Nele,MUA.nod);
dhdBnod=reshape(dhdB(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F.q)
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    qnod=[];
end

if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    muknod=[];
end

if ~isempty(F.V0)
    V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    V0nod=[];
end

T=zeros(MUA.Nele,MUA.nod);

CtrlVar.BasalDrag.CalculateDerivatives=true;   % dtaubxdh and dtaubydh are needed below

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);

    Dx=reshape(Deriv(:,1,:),MUA.Nele,MUA.nod);
    Dy=reshape(Deriv(:,2,:),MUA.Nele,MUA.nod);

    % values at the integration point

    hint=hnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;

    uint=unod*fun;
    vint=vnod*fun;

    mint=mnod*fun;
    nint=nnod*fun;

    % C and AGlen are coefficients here, B is the differentiation variable, so no chain-rule factor is needed
    Cint=SmoothFloor(Cnod*fun,CtrlVar.Cmin,CtrlVar.CminWidth);
    AGlenInt=SmoothFloor(AGlennod*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

    if ~isempty(qnod)    ; qint=qnod*fun      ; else ; qint=[]    ; end
    if ~isempty(muknod)  ; mukint=muknod*fun  ; else ; mukint=[]  ; end
    if ~isempty(V0nod)   ; V0int=V0nod*fun    ; else ; V0int=[]   ; end

    Psi_x_int=Psi_xnod*fun;
    Psi_y_int=Psi_ynod*fun;

    % Masks. Because the Heaviside is symmetrical, the complementary mask is 1-Heint and the complementary delta is
    % deltaint, so neither needs to be formed separately.
    hfint=F.rhow*Hint./rhoint;
    Heint    = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    deltaint = DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);

    HeHint    = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0);
    deltaHint = DiracDelta(CtrlVar.kH,Hint,CtrlVar.Hh0);
    Hposint   = HeHint.*Hint;

    dint = (1-Heint).*rhoint.*hint/F.rhow + Heint.*Hposint ;

    % gradients at the integration point
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

    if isempty(dtaubxdh) || isempty(dtaubydh)
        error("dIdBqGeneral:NoThicknessDerivatives",...
            "BasalDrag.m did not return dtaubxdh/dtaubydh for the sliding law %s.\n",CtrlVar.SlidingLaw)
    end

    detJw=detJ*MUA.weights(Iint);

    % Building blocks of the derivative of the draft, d.
    % Theta is delta multiplied by the difference between the two candidate drafts, and therefore vanishes at flotation.
    Theta  = deltaint.*( Hposint - rhoint.*hint/F.rhow ) ;   % dd/d(Dh)
    dddhF  = (1-Heint).*rhoint/F.rhow ;                      % dd/dh at fixed Dh
    dddB_H = Heint.*( deltaHint.*Hint + HeHint ) ;            % -dd/dB through H=S-B

    % Contractions with the adjoint fields, all independent of the free node index
    Visc = etaint.*( (4*exx+2*eyy).*dlxdx + 2*exy.*dlxdy + (4*eyy+2*exx).*dlydy + 2*exy.*dlydx ) ;
    Drag = dtaubxdh.*Psi_x_int + dtaubydh.*Psi_y_int ;
    Grav = dbdx.*Psi_x_int + dbdy.*Psi_y_int ;
    Pres = dlxdx + dlydy ;
    Buoy = rhoint.*hint - F.rhow*dint ;

    for Inod=1:MUA.nod

        hp = dhdBnod(:,Inod) ;              % (dh/dB)_j
        bp = dbdBnod(:,Inod) ;              % (db/dB)_j

        mu = hp + F.rhow./rhoint ;          % (d(Dh)/dB)_j , the flotation measure

        dd = Theta.*mu + dddhF.*hp - dddB_H ;   % (dd/dB)_j , as a coefficient of phi_j

        % everything multiplying phi_j
        A =   hp.*Visc ...
            + mu.*Drag ...
            + ca*F.g*( rhoint.*hp - F.rhow*dd ).*Grav ...
            - F.g*sa*rhoint.*hp.*Psi_x_int ...
            - ca*F.g*( rhoint.*hint.*hp - F.rhow*dint.*dd ).*Pres ;

        % everything multiplying the gradient of phi_j, from d(db/dx)/dB_j = (db/dB)_j dphi_j/dx
        Agrad = ca*F.g*Buoy.*bp ;

        T(:,Inod)=T(:,Inod) ...
            + ( A.*fun(Inod) + Agrad.*( Dx(:,Inod).*Psi_x_int + Dy(:,Inod).*Psi_y_int ) ).*detJw ;

    end
end

dIdB=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdB=dIdB+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

% No conversion to log space: B is inverted for directly, and can be negative.

%% test the gradient

if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdB);
end

end




function FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dIdB)  %#ok<INUSD>

nTests=150; 
iNodeVector=randi(MUA.Nnodes,1,nTests);
dIdB_vector=nan(nTests,1);
dIdB_FD_vector=nan(nTests,1);

B0=F.B;

% The grounding-line transition has a width of about 1/kH, so the step must be small compared with that.
dB=1e-3/max(CtrlVar.kH,1e-10) ;

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=true;
CtrlVar.Calculate.Geometry="bs-FROM-hBS";

if ~isfield(CtrlVar,"MapOldToNew") || ~isfield(CtrlVar.MapOldToNew,"Test")
    CtrlVar.MapOldToNew.Test=false;
end


for iTest=1:numel(iNodeVector)
    iNode=iNodeVector(iTest);

    F.B=B0 ; F.B(iNode)=B0(iNode)-dB ;
    [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow);
    Ruv_minus = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
    b_minus = Ruv_minus.' * [Psi_x; Psi_y];

    F.B=B0 ; F.B(iNode)=B0(iNode)+dB ;
    [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow);
    Ruv_plus = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
    b_plus = Ruv_plus.' * [Psi_x; Psi_y];

    F.B=B0 ;

    dIdB_FD = (b_plus - b_minus)/(2*dB);

    dIdB_FD_vector(iTest)= dIdB_FD;
    dIdB_vector(iTest)=dIdB(iNode);

    Diff=norm(dIdB(iNode) - dIdB_FD)/(abs(dIdB(iNode))+eps);
    if nTests<2
        fprintf("dIdBqGeneral: normalized norm of difference between Direct-Adjoint and FD for node %i is %g \n",iNode,Diff)
    end
end

if nTests>2

    figdIDB=FindOrCreateFigure("Test:dIdB") ; clf(figdIDB)


    fig_dIdBTest=FindOrCreateFigure("Test dIdB") ; clf(fig_dIdBTest)
    plot(dIdB_FD_vector,dIdB_vector,"or") 
    axis equal
    hold on ;
    plot([min(dIdB_vector) max(dIdB_vector)],[min(dIdB_vector) max(dIdB_vector)],"--k")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off
    xlabel("$\langle  \delta_{B_i} F^x \phi_i | \Psi_x \rangle  $",Interpreter="latex")  ;
    ylabel("Finite differences",Interpreter="latex")
    title("$ \langle  \delta_{B_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{B_i} F^y \phi_i| \Psi_y \rangle $",Interpreter="latex")
    subtitle(sprintf("Normalized diff %g",Diff),Interpreter="latex")

    

end

end
