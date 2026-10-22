function KJBq=JhdotBq(CtrlVar,MUA,F,Meas)

narginchk(4,4)
nargoutchk(1,1)

%% Explicit mixed second derivative of the dh/dt misfit term with respect to B and q=(u,v)
%
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left ( \frac{\dot{h} - \dot{h}_{Meas}}{\dot{h}_{err}} \right )^2 dx \, dy ,  \qquad \dot{h} = a - \frac{1}{\rho} \nabla \cdot (\rho h \mathbf{v}) $$
%
% The model hdot depends on the velocities and on the thickness, and (with the upper surface s held fixed) the thickness is a function of B,
% h=h(B). J_hdot is therefore an explicit function of both q and B, and the mixed second derivative J^{pq}, with p=B, is not zero.
%
% With H'' = diag(dh/dB) (from dGeometrydB.m),
%
% $$ J^{Bq} = H'' \, \left [ J_{hu} \;\; J_{hv} \right ] , \qquad (J_{hu})_{ji} = \frac{\partial^2 J_{\dot{h}}}{\partial u_i \partial h_j} $$
%
% $$ \frac{\partial^2 J_{\dot{h}}}{\partial u_i \partial h_j} = \frac{1}{\mathcal{A}} \int \! \int \dot{h}_{err}^{-2} \left [ G^x_i \, S_j - (\dot{h}-\dot{h}_{Meas}) \, T^x_{ij} \right ] dx \, dy $$
%
% with
%
% $$ G^x_i = \left ( \partial_x h + h \frac{\partial_x \rho}{\rho} \right ) \phi_i + h \, \partial_x \phi_i  $$
%
% $$ S_j = \left ( \partial_x \phi_j + \phi_j \frac{\partial_x \rho}{\rho} \right ) u + \phi_j \partial_x u + \left ( \partial_y \phi_j + \phi_j \frac{\partial_y \rho}{\rho} \right ) v + \phi_j \partial_y v $$
%
% $$ T^x_{ij} = \left ( \partial_x \phi_j + \phi_j \frac{\partial_x \rho}{\rho} \right ) \phi_i + \phi_j \, \partial_x \phi_i $$
%
% and the v-derivative with x replaced by y. The matrix returned is n x 2n, with the row index the bedrock node and the column index the
% velocity node, ordered (u,v), i.e. in the same layout as the matrices returned by FBuv.m.
%
% Meas.dhdtCov is assumed to be diagonal, as in Jq.m, Jqq.m and Misfit.m.
%
% see also: Hess_qp.m, JhdotBB.m, FBuv.m, EvaluateJhdotAndDerivatives.m, dGeometrydB.m
%%

ndim=2 ; Nele=MUA.Nele ; nod=MUA.nod ; Nnodes=MUA.Nnodes ; Area=MUA.Area ;

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA) ;
end

anod=reshape(F.as(MUA.connectivity,1),Nele,nod)+reshape(F.ab(MUA.connectivity,1),Nele,nod) ;
hnod=reshape(F.h(MUA.connectivity,1),Nele,nod) ;
unod=reshape(F.ub(MUA.connectivity,1),Nele,nod) ;
vnod=reshape(F.vb(MUA.connectivity,1),Nele,nod) ;
rhonod=reshape(F.rho(MUA.connectivity,1),Nele,nod) ;
dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),Nele,nod) ;

dhdtErr=full(sqrt(spdiags(Meas.dhdtCov))) ;
dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),Nele,nod) ;

Ju=zeros(Nele,nod,nod) ;      % stored as [Nele , j (thickness node) , i (velocity node)]
Jv=zeros(Nele,nod,nod) ;

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint) ;
    Deriv=MUA.Deriv(:,:,:,Iint) ;
    Dx=reshape(Deriv(:,1,:),Nele,nod) ;
    Dy=reshape(Deriv(:,2,:),Nele,nod) ;

    aint=anod*fun ; hint=hnod*fun ; uint=unod*fun ; vint=vnod*fun ; rhoint=rhonod*fun ;
    hdotMeasint=dhdtMeasnod*fun ; errint=dhdtErrnod*fun ;

    dhdx=sum(Dx.*hnod,2) ; dhdy=sum(Dy.*hnod,2) ;
    dudx=sum(Dx.*unod,2) ; dvdy=sum(Dy.*vnod,2) ;
    drhodx=sum(Dx.*rhonod,2) ; drhody=sum(Dy.*rhonod,2) ;

    % the model dh/dt at the integration point, as in EvaluateJhdotAndDerivatives.m
    hdot=aint-(rhoint.*dhdx.*uint+rhoint.*hint.*dudx+drhodx.*hint.*uint+rhoint.*dhdy.*vint+rhoint.*hint.*dvdy+drhody.*hint.*vint)./rhoint ;
    res=hdot-hdotMeasint ;

    W=detJ.*MUA.weights(Iint)./(Area.*errint.^2) ;

    kx=dhdx+hint.*drhodx./rhoint ; ky=dhdy+hint.*drhody./rhoint ;
    funR=fun.' ;
    Gx=kx.*funR+hint.*Dx ;  Gy=ky.*funR+hint.*Dy ;                 % Nele x nod, index i
    Px=Dx+(drhodx./rhoint).*funR ;  Py=Dy+(drhody./rhoint).*funR ;  % Nele x nod, index j
    S=Px.*uint+funR.*dudx+Py.*vint+funR.*dvdy ;                     % Nele x nod, index j

    for Inod=1:nod
        Tx=Px.*fun(Inod)+funR.*Dx(:,Inod) ;    % Nele x nod (index j)
        Ty=Py.*fun(Inod)+funR.*Dy(:,Inod) ;
        Ju(:,:,Inod)=Ju(:,:,Inod)+W.*(Gx(:,Inod).*S-res.*Tx) ;
        Jv(:,:,Inod)=Jv(:,:,Inod)+W.*(Gy(:,Inod).*S-res.*Ty) ;
    end

end

nnzTot=nod*nod*Nele ; Iind=zeros(nnzTot,1) ; Jind=zeros(nnzTot,1) ; Xu=zeros(nnzTot,1) ; Xv=zeros(nnzTot,1) ; istak=0 ;
for Inod=1:nod
    for Jnod=1:nod
        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Jnod) ;     % row: thickness (bedrock) node j
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Inod) ;     % column: velocity node i
        Xu(istak+1:istak+Nele)=Ju(:,Jnod,Inod) ;
        Xv(istak+1:istak+Nele)=Jv(:,Jnod,Inod) ;
        istak=istak+Nele ;
    end
end
Jhu=sparse(Iind,Jind,Xu,Nnodes,Nnodes) ;
Jhv=sparse(Iind,Jind,Xv,Nnodes,Nnodes) ;

[~,dhdB]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow) ;
Dh=spdiags(dhdB(:),0,Nnodes,Nnodes) ;

KJBq=[Dh*Jhu , Dh*Jhv] ;

end
