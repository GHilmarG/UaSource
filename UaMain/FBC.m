function KFBC=FBC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)  %#ok<INUSD>

narginchk(7,7)
nargoutchk(1,1)

%% Builds the Hessian cross term
%
% $$ \mathcal{F}^{BC} = \langle \Psi , \delta^2_{BC} \mathcal{F} \rangle $$
%
% i.e. the second-order mixed derivative of the forward model with respect to the bedrock, B, and the slipperiness, C, contracted with the
% adjoint variables. The matrix returned is n x n, with the row index the B node and the column index the C node.
%
%% Structure of the calculation
%
% C enters the forward model only through the basal traction, $t_b = \mathcal{G} \, \beta^2(C) \, \mathbf{u}$, where $\mathcal{G}$ is the
% grounding (flotation) mask, $\mathcal{G} = \mathcal{H}(\Delta h)$ with $\Delta h = h - h_f$. For the Weertman sliding law $\beta^2$ does not depend on
% B (nor on h), and B enters the traction only through $\mathcal{G}$.
%
% As in the forward model (uvMatrixAssemblySSTREAM.m), and in dIdBqGeneral.m, FBuv.m and FBB.m, the flotation thickness is evaluated at the
% integration point, $h_f = \rho_o (S-B)/\rho$ with $S$, $B$ and $\rho$ interpolated to the integration point. Then
%
% $$ \frac{\partial \Delta h}{\partial B_j} = \left ( h'_j + \frac{\rho_o}{\rho} \right ) \phi_j $$
%
% where $h'$ is from dGeometrydB.m and $\rho$ is the density at the integration point. With $\partial \mathcal{G} / \partial \Delta h = \delta(\Delta h)$,
%
% $$ \mathcal{F}^{BC}_{jk} = \int \! \int \delta(\Delta h) \left ( h'_j + \frac{\rho_o}{\rho} \right ) \frac{\partial \beta^2}{\partial C} \frac{d C_{eff}}{d C}
%    \, ( u \Psi_x + v \Psi_y ) \; \phi_j \phi_k \; dx \, dy $$
%
% which is assembled as diag($h'$) times a weighted mass matrix, plus a second weighted mass matrix with the weight multiplied by $\rho_o/\rho$.
% The kernel is that of dIdCq.m with the grounding mask replaced by its derivative. If the inversion is for log10(C) the columns are
% multiplied by $C \ln 10$.
%
% Note: dIdCq.m and FCC.m evaluate $h_f$ from the nodal values, which is only the same as the integration-point evaluation when the density
% is uniform.
%
% Only implemented for the Weertman sliding law, as FCuv.m and Fqq.m.
%
% see also: FCC.m, FBB.m, FBuv.m, FAB.m, dIdCq.m, dGeometrydB.m, Fpp.m
%%

KFBC=[] ;

[~,isB,isC]=isABC(CtrlVar) ;

if ~(isB && isC)
    % Returning an empty matrix allows Fpp.m to treat inactive fields uniformly
    return
end

if CtrlVar.IncludeMelangeModelPhysics
    error("FBC:MelangeNotImplemented","FBC is not implemented for melange model physics.")
end

if CtrlVar.Hh0~=0
    error("FBC:NonZeroHh0","A symmetrical Heaviside function is assumed, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

if ~ismember(string(CtrlVar.SlidingLaw),["W","Weertman"])
    error("FBC:SlidingLawNotImplemented","FBC is implemented for the Weertman sliding law only, but CtrlVar.SlidingLaw=""%s""",string(CtrlVar.SlidingLaw))
end

ndim=2 ; nNodes=MUA.Nnodes ; Nele=MUA.Nele ; nod=MUA.nod ; nip=MUA.nip ;
C0=CtrlVar.Czero ; u0=CtrlVar.SpeedZero ;

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA) ;
end

h_node=reshape(F.h(MUA.connectivity,1),Nele,nod) ;
B_node=reshape(F.B(MUA.connectivity,1),Nele,nod) ;
S_node=reshape(F.S(MUA.connectivity,1),Nele,nod) ;
rho_node=reshape(F.rho(MUA.connectivity,1),Nele,nod) ;

C_node=reshape(F.C(MUA.connectivity,1),Nele,nod) ;
m_node=reshape(F.m(MUA.connectivity,1),Nele,nod) ;
u_node=reshape(F.ub(MUA.connectivity,1),Nele,nod) ;
v_node=reshape(F.vb(MUA.connectivity,1),Nele,nod) ;
Psi_x_node=reshape(Psi_x(MUA.connectivity,1),Nele,nod) ;
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),Nele,nod) ;

W1=zeros(Nele,nip) ;         % integration-point weights, multiplying h'_j
W2=zeros(Nele,nip) ;         % integration-point weights, multiplying rho_o/rho
P=zeros(nip,nod*nod) ;       % outer products of the shape functions

for Iint=1:nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint) ;

    h=h_node*fun ; Bint=B_node*fun ; Sint=S_node*fun ; rhoint=rho_node*fun ;
    hfint=F.rhow*(Sint-Bint)./rhoint ;                   % as in uvMatrixAssemblySSTREAM.m
    u=u_node*fun ; v=v_node*fun ;
    u_e=sqrt(u.*u+v.*v+u0^2) ;

    [C,dCeffdC]=SmoothFloor(C_node*fun,CtrlVar.Cmin,CtrlVar.CminWidth) ;      % smooth clipping, as in FCC.m and dIdCq.m
    m=m_node*fun ;
    Psix=Psi_x_node*fun ; Psiy=Psi_y_node*fun ;

    delta=DiracDelta(CtrlVar.kH,h-hfint,CtrlVar.Hh0) ;     % d G / d (h-hf)

    CC=C+C0 ;
    tBeta=(u_e./CC).^(1./m) ;
    dbeta2dC=-(1./m).*tBeta./(CC.*u_e) ;               % d beta^2 / d C  (Weertman)

    Kernel=delta.*dbeta2dC.*dCeffdC.*(Psix.*u+Psiy.*v) ;

    detJw=detJ*MUA.weights(Iint) ;
    W1(:,Iint)=Kernel.*detJw ;
    W2(:,Iint)=Kernel.*detJw.*(F.rhow./rhoint) ;
    P(Iint,:)=reshape(fun*fun.',1,nod*nod) ;

end

H1=reshape(W1*P,Nele,nod,nod) ;
H2=reshape(W2*P,Nele,nod,nod) ;

% collect all values and indices into vectors and only make one call to the sparse function
Iind=zeros(nod*nod*Nele,1,'uint32') ; Jind=zeros(nod*nod*Nele,1,'uint32') ; X1=zeros(nod*nod*Nele,1) ; X2=zeros(nod*nod*Nele,1) ; istak=0 ;
for Inode=1:nod
    for Jnode=1:nod
        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inode) ;
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Jnode) ;
        X1(istak+1:istak+Nele)=H1(:,Inode,Jnode) ;
        X2(istak+1:istak+Nele)=H2(:,Inode,Jnode) ;
        istak=istak+Nele ;
    end
end
M1=sparse(Iind,Jind,X1,nNodes,nNodes) ;
M2=sparse(Iind,Jind,X2,nNodes,nNodes) ;

%% the B dependence enters through the flotation measure, d(h-hf)/dB_j = (dh/dB_j + rho_o/rho) phi_j
[~,dhdB]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow) ;
KFBC=spdiags(dhdB(:),0,nNodes,nNodes)*M1+M2 ;

%% from C to log10(C) (B is not transformed)
if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    KFBC=KFBC*spdiags(F.C(:)*log(10),0,nNodes,nNodes) ;
end

end
