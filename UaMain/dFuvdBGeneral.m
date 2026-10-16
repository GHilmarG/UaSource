


function K=dFuvdBGeneral(CtrlVar,MUA,F)

%% Calculates the sensitivity matrix
%
% $$ \partial F^{uv} / \partial B $$
%
% If $n$ is the number of nodes, the matrix returned is $2 n \times n$
%
% $$\left[\begin{array}{cccc}
% \partial F^x_1 /\partial B_1  & \partial F^x_1 /\partial B_2  & \ldots & \partial F^x_1 /\partial B_n  \\
% \partial F^x_2 /\partial B_1  & \partial F^x_2 /\partial B_2  & \ldots & \partial F^x_2 /\partial B_n  \\
%              .               &              .               &  .     &    .                          \\
% \partial F^y_1 /\partial B_1  & \partial F^y_1 /\partial B_2  & \ldots & \partial F^y_1 /\partial B_n  \\
%              .               &              .               &  .     &    .                          \\
% \end{array}\right] $$
%
% This is the general version of
%
%   dFuvdB.m
%
% Unlike dFuvdB.m, no assumption is made about the ice being grounded. The floating mask is itself a function of B and
% its variation with B is accounted for.
%
%% Relationship to dIdBqGeneral.m
%
% This matrix is the un-contracted form of the quantity returned by dIdBqGeneral.m, i.e.
%
% $$ \mathrm{dIdBqGeneral} = \left ( \frac{\partial F^{uv}}{\partial B} \right )^{T} \left [ \Psi_x ; \Psi_y \right ] $$
%
% Every term below is taken directly from dIdBqGeneral.m, the only difference being that the free index j is retained
% rather than contracted against the adjoint fields. Where dIdBqGeneral.m multiplies by Psi_x_int, this function
% multiplies by fun(Inod), and where dIdBqGeneral.m multiplies by dlxdx, this function multiplies by Dx(:,Inod).
%
% The above identity is the recommended test of this function: see TestdFuvdBGeneral.m
%
%% Assumptions
%
% As in dIdBqGeneral.m:
%
% # The upper ice surface, s, the ocean surface, S, and the densities are held fixed, and b and h are obtained from s,
%   S and B using the closure solved by Calc_bh_From_sBS.m. F.b and F.h must be the converged output of that closure
%   for the current F.B.
% # The Heaviside function is symmetrical, i.e. CtrlVar.Hh0=0.
% # All sliding laws are supported. Melange model physics is not.
%
%% The terms
%
% Writing the nodal geometry derivatives as
%
% $$ b'_j = (\partial b/\partial B)_j , \qquad h'_j = -b'_j , \qquad \mu_j = h'_j + \rho_o/\rho $$
%
% and the derivative of the draft as
%
% $$ d'_j = \Theta \mu_j + (1-\mathcal{G}) \frac{\rho}{\rho_o} h'_j - \mathcal{G} ( \delta(H) H + \mathcal{H}(H) ) ,
%    \qquad \Theta = \delta(\Delta h) \left ( H^{+} - \frac{\rho h}{\rho_o} \right ) $$
%
% the x-component is
%
% $$ \frac{\partial F^x_i}{\partial B_j} =
%  \langle h'_j \eta ( 4 \partial_x u + 2 \partial_y v) \phi_j | \partial_x \phi_i \rangle
% + \langle h'_j \eta ( \partial_y u + \partial_x v) \phi_j | \partial_y \phi_i \rangle
% + \langle \mu_j \, \partial_h t_{bx} \, \phi_j | \phi_i \rangle $$
%
% $$ + \, g \cos\alpha \langle ( \rho h'_j - \rho_o d'_j ) \phi_j \, \partial_x b
%      + (\rho h - \rho_o d ) \, b'_j \, \partial_x \phi_j  | \phi_i \rangle $$
%
% $$ - \, \rho g \sin\alpha \langle h'_j \phi_j | \phi_i \rangle
%    - \, g \cos\alpha \langle ( \rho h h'_j - \rho_o d d'_j ) \phi_j | \partial_x \phi_i \rangle $$
%
% and correspondingly for $F^y_i$, with x and y interchanged in the differentiation slots and without the
% $\sin \alpha$ term.
%
%  see also: dIdBqGeneral.m, dGeometrydB.m, dFuvdB.m, duvdBGeneral.m
%
%%

narginchk(3,3)
nargoutchk(1,1)

if CtrlVar.IncludeMelangeModelPhysics
    error("dFuvdBGeneral:MelangeNotImplemented","Sensitivities to B with melange model physics are not implemented.")
end

if CtrlVar.Hh0~=0
    error("dFuvdBGeneral:NonZeroHh0",...
        "This implementation assumes a symmetrical Heaviside function, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

ndim=2;
ca=cos(F.alpha) ; sa=sin(F.alpha) ;
nNodes=MUA.Nnodes ;

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

dbdBnod=reshape(dbdB(MUA.connectivity,1),MUA.Nele,MUA.nod);
dhdBnod=reshape(dhdB(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F.q)   ; qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod)     ; else ; qnod=[]   ; end
if ~isempty(F.muk) ; muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod) ; else ; muknod=[] ; end
if ~isempty(F.V0)  ; V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod)   ; else ; V0nod=[]  ; end

dFxdB=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFydB=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end

CtrlVar.BasalDrag.CalculateDerivatives=true;   % dtaubxdh and dtaubydh are needed below

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

    % C and AGlen are coefficients here, B is the differentiation variable, so no chain-rule factor is needed
    Cint=SmoothFloor(Cnod*fun,CtrlVar.Cmin,CtrlVar.CminWidth);
    AGlenInt=SmoothFloor(AGlennod*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

    if ~isempty(qnod)   ; qint=qnod*fun     ; else ; qint=[]   ; end
    if ~isempty(muknod) ; mukint=muknod*fun ; else ; mukint=[] ; end
    if ~isempty(V0nod)  ; V0int=V0nod*fun   ; else ; V0int=[]  ; end

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

    etaint=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);

    [~,~,~,~,~,~,dtaubxdh,dtaubydh] = ...
        BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,...
        [],[],[],[],[],[],[],[],qint,F.g,mukint,V0int);

    if isempty(dtaubxdh) || isempty(dtaubydh)
        error("dFuvdBGeneral:NoThicknessDerivatives",...
            "BasalDrag.m did not return dtaubxdh/dtaubydh for the sliding law %s.\n",CtrlVar.SlidingLaw)
    end

    detJw=detJ*MUA.weights(Iint);

    % Building blocks of the derivative of the draft, d.
    % Theta is delta multiplied by the difference between the two candidate drafts, and therefore vanishes at flotation.
    Theta  = deltaint.*( Hposint - rhoint.*hint/F.rhow ) ;   % dd/d(Dh)
    dddhF  = (1-Heint).*rhoint/F.rhow ;                      % dd/dh at fixed Dh
    dddB_H = Heint.*( deltaHint.*Hint + HeHint ) ;           % -dd/dB through H=S-B

    Buoy = rhoint.*hint - F.rhow*dint ;

    for Jnod=1:MUA.nod     % the parameter index, i.e. the node at which B is perturbed

        hp = dhdBnod(:,Jnod) ;              % (dh/dB)_j
        bp = dbdBnod(:,Jnod) ;              % (db/dB)_j
        mu = hp + F.rhow./rhoint ;          % (d(Dh)/dB)_j , the flotation measure
        dd = Theta.*mu + dddhF.*hp - dddB_H ;   % (dd/dB)_j , as a coefficient of phi_j

        % terms which, in the x equation, multiply phi_i
        Ax = ( mu.*dtaubxdh ...
             + ca*F.g*( rhoint.*hp - F.rhow*dd ).*dbdx ...
             - F.g*sa*rhoint.*hp ).*fun(Jnod) ...
             + ca*F.g*Buoy.*bp.*Dx(:,Jnod) ;

        Ay = ( mu.*dtaubydh ...
             + ca*F.g*( rhoint.*hp - F.rhow*dd ).*dbdy ).*fun(Jnod) ...
             + ca*F.g*Buoy.*bp.*Dy(:,Jnod) ;

        % terms which multiply the gradient of phi_i
        Pref = -ca*F.g*( rhoint.*hint.*hp - F.rhow*dint.*dd ).*fun(Jnod) ;   % pressure term
        Vxx  = hp.*etaint.*(4*exx+2*eyy).*fun(Jnod) ;    % viscous, x eq, d_x phi_i
        Vxy  = hp.*etaint.*(2*exy).*fun(Jnod) ;          % viscous, x eq, d_y phi_i
        Vyy  = hp.*etaint.*(4*eyy+2*exx).*fun(Jnod) ;    % viscous, y eq, d_y phi_i
        Vyx  = Vxy ;                                     % viscous, y eq, d_x phi_i

        for Inod=1:MUA.nod   % the test-function index

            dFxdB(:,Inod,Jnod)=dFxdB(:,Inod,Jnod) ...
                + ( Ax.*fun(Inod) + (Vxx+Pref).*Dx(:,Inod) + Vxy.*Dy(:,Inod) ).*detJw ;

            dFydB(:,Inod,Jnod)=dFydB(:,Inod,Jnod) ...
                + ( Ay.*fun(Inod) + (Vyy+Pref).*Dy(:,Inod) + Vyx.*Dx(:,Inod) ).*detJw ;

        end
    end
end

%% assemble into a 2n x n sparse matrix

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1);

istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=dFxdB(:,Inod,Jnod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+nNodes;
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=dFydB(:,Inod,Jnod);
        istak=istak+MUA.Nele;

    end
end

K=sparse(Iind,Jind,Xval,2*nNodes,nNodes);

end
