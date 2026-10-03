function KFAB=FAB(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)  %#ok<INUSD>

narginchk(7,7)
nargoutchk(1,1)

%% Builds the Hessian cross term
%
% $$ \mathcal{F}^{AB} = \langle \Psi , \delta^2_{AB} \mathcal{F} \rangle $$
%
% i.e. the second-order mixed derivative of the forward model with respect to the rate factor, A, and the bedrock, B, contracted with the
% adjoint variables. The matrix returned is n x n, with the row index the A node and the column index the B node.
%
%% Structure of the calculation
%
% As described in FBB.m, the contracted residual density contains the viscous part $h \, V$, with
%
% $$ V = \eta \Big[ (4 \partial_x u + 2\partial_y v)\partial_x \Psi_x + (\partial_y u + \partial_x v)\partial_y \Psi_x
%        + (4 \partial_y v + 2\partial_x u)\partial_y \Psi_y + (\partial_y u + \partial_x v)\partial_x \Psi_y \Big] $$
%
% A enters the forward model only through the viscosity $\eta$, and B enters this term only through the thickness, $h=h(B)$,
% $\partial h / \partial B_j = h'_j \phi_j$ with $h'$ from dGeometrydB.m. The sliding law, the driving-stress terms and the
% grounding mask do not depend on A, and $\eta$ does not depend on B. Hence
%
% $$ \mathcal{F}^{AB}_{kj} = h'_j \int \! \int \frac{\partial \eta}{\partial A} \frac{d A_{eff}}{d A} \, T \; \phi_k \phi_j \; dx \, dy , \qquad
%    T = (4 \partial_x u + 2\partial_y v)\partial_x \Psi_x + 2 e_{xy} \partial_y \Psi_x + (4 \partial_y v + 2\partial_x u)\partial_y \Psi_y + 2 e_{xy} \partial_x \Psi_y $$
%
% which is a mass matrix weighted at the integration points, times diag(h') on the B index. This is the same kernel as in dIdAq.m, with h
% replaced by $h'_j \phi_j$. If the inversion is for log10(A), the rows are multiplied by $A \ln 10$ (B is a linear variable, and A enters
% locally, so no further terms arise).
%
% see also: FAA.m, FBB.m, FBuv.m, FBC.m, dIdAq.m, dGeometrydB.m, Fpp.m
%%

KFAB=[] ;

[isA,isB,~]=isABC(CtrlVar) ;

if ~(isA && isB)
    % Returning an empty matrix allows Fpp.m to treat inactive fields uniformly
    return
end

if CtrlVar.IncludeMelangeModelPhysics
    error("FAB:MelangeNotImplemented","FAB is not implemented for melange model physics.")
end

if CtrlVar.Hh0~=0
    error("FAB:NonZeroHh0","A symmetrical Heaviside function is assumed, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

ndim=2 ; nNodes=MUA.Nnodes ; Nele=MUA.Nele ; nod=MUA.nod ; nip=MUA.nip ;

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA) ;
end

A_node=reshape(F.AGlen(MUA.connectivity,1),Nele,nod) ;
n_node=reshape(F.n(MUA.connectivity,1),Nele,nod) ;
u_node=reshape(F.ub(MUA.connectivity,1),Nele,nod) ;
v_node=reshape(F.vb(MUA.connectivity,1),Nele,nod) ;
Psi_x_node=reshape(Psi_x(MUA.connectivity,1),Nele,nod) ;
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),Nele,nod) ;

W=zeros(Nele,nip) ;          % integration-point weights
P=zeros(nip,nod*nod) ;       % outer products of the shape functions

CtrlVar.EffectiveViscosity.CalculateDerivatives=false ;

for Iint=1:nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint) ;
    detJ=MUA.DetJ(:,Iint) ;
    Dx=reshape(Deriv(:,1,:),Nele,nod) ;
    Dy=reshape(Deriv(:,2,:),Nele,nod) ;

    n=n_node*fun ;
    [A,dAeffdA]=SmoothFloor(A_node*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth) ;

    exx=sum(Dx.*u_node,2) ; eyy=sum(Dy.*v_node,2) ; exy=0.5*(sum(Dy.*u_node,2)+sum(Dx.*v_node,2)) ;

    dPsi_x_dx=sum(Dx.*Psi_x_node,2) ; dPsi_x_dy=sum(Dy.*Psi_x_node,2) ;
    dPsi_y_dx=sum(Dx.*Psi_y_node,2) ; dPsi_y_dy=sum(Dy.*Psi_y_node,2) ;

    [~,~,~,detadA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy) ;
    detadA=detadA.*dAeffdA ;

    Temp=(4*exx+2*eyy).*dPsi_x_dx + 2*exy.*dPsi_x_dy + (4*eyy+2*exx).*dPsi_y_dy + 2*exy.*dPsi_y_dx ;

    detJw=detJ*MUA.weights(Iint) ;
    W(:,Iint)=detadA.*Temp.*detJw ;
    P(Iint,:)=reshape(fun*fun.',1,nod*nod) ;

end

H=reshape(W*P,Nele,nod,nod) ;

% collect all values and indices into vectors and only make one call to the sparse function
Iind=zeros(nod*nod*Nele,1,'uint32') ; Jind=zeros(nod*nod*Nele,1,'uint32') ; Xval=zeros(nod*nod*Nele,1) ; istak=0 ;
for Inode=1:nod
    for Jnode=1:nod
        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inode) ;
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Jnode) ;
        Xval(istak+1:istak+Nele)=H(:,Inode,Jnode) ;
        istak=istak+Nele ;
    end
end
M=sparse(Iind,Jind,Xval,nNodes,nNodes) ;

%% the B dependence enters through the nodal thickness derivative dh/dB
[~,dhdB]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow) ;
KFAB=M*spdiags(dhdB(:),0,nNodes,nNodes) ;

%% from A to log10(A) (B is not transformed)
if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    KFAB=spdiags(F.AGlen(:)*log(10),0,nNodes,nNodes)*KFAB ;
end

end
