function KhdotBB=JhdotBB(CtrlVar,MUA,F,Meas)

narginchk(4,4)
nargoutchk(1,1)

%% Explicit second derivative, with respect to B, of the dh/dt misfit term (velocities held fixed)
%
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left ( \frac{\dot{h} - \dot{h}_{Meas}}{\dot{h}_{err}} \right )^2 dx \, dy ,  \qquad \dot{h} = a - \frac{1}{\rho} \nabla \cdot (\rho h \mathbf{v}) $$
%
% The model hdot depends on the thickness, and when inverting for B (with the upper surface s held fixed) the thickness is a function of B, h=h(B),
% through the geometrical closure. J_hdot is therefore an explicit function of B, i.e. it depends on B other than through u and v.
%
% With H'' = diag(dh/dB) and h2 = d^2h/dB^2 (from dGeometrydB.m), the explicit second derivative with respect to B is
%
% $$ J^{BB}_{\dot{h}} = H'' \, J_{hh} \, H'' + \mathrm{diag} ( J_h \, h2 ) $$
%
% where J_h is the derivative of J_hdot with respect to the nodal thickness (the fourth output of EvaluateJhdotAndDerivatives.m), and
%
% $$ (J_{hh})_{ij} = \frac{1}{\mathcal{A}} \int \! \int \dot{h}_{err}^{-2} \, S_i \, S_j \; dx \, dy $$
%
% $$ S_j = \left ( \partial_x \phi_j + \phi_j \frac{\partial_x \rho}{\rho} \right ) u + \phi_j \partial_x u + \left ( \partial_y \phi_j + \phi_j \frac{\partial_y \rho}{\rho} \right ) v + \phi_j \partial_y v $$
%
% There is no term proportional to the residual in J_hh because hdot is linear in h.
%
% Meas.dhdtCov is assumed to be diagonal, as in Jq.m, Jqq.m and Misfit.m.
%
% see also: Jpp.m, EvaluateJhdotAndDerivatives.m, dGeometrydB.m, Jqq.m
%%

ndim=2 ; Nele=MUA.Nele ; nod=MUA.nod ; Nnodes=MUA.Nnodes ; Area=MUA.Area ;

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA) ;
end

unod=reshape(F.ub(MUA.connectivity,1),Nele,nod) ;
vnod=reshape(F.vb(MUA.connectivity,1),Nele,nod) ;
rhonod=reshape(F.rho(MUA.connectivity,1),Nele,nod) ;

dhdtErr=full(sqrt(spdiags(Meas.dhdtCov))) ;
dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),Nele,nod) ;

Jhh=zeros(Nele,nod,nod) ;      % stored as [Nele , k , i]

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint) ;
    Deriv=MUA.Deriv(:,:,:,Iint) ;
    Dx=reshape(Deriv(:,1,:),Nele,nod) ;
    Dy=reshape(Deriv(:,2,:),Nele,nod) ;

    uint=unod*fun ; vint=vnod*fun ; rhoint=rhonod*fun ; errint=dhdtErrnod*fun ;

    dudx=sum(Dx.*unod,2) ;  dvdy=sum(Dy.*vnod,2) ;
    drhodx=sum(Dx.*rhonod,2) ;  drhody=sum(Dy.*rhonod,2) ;

    W=detJ.*MUA.weights(Iint)./(Area.*errint.^2) ;     % quadrature weight folded in

    % S(:,j) is the variation of hdot with respect to the thickness at node j (up to a sign), Nele x nod
    funR=fun.' ;
    S=(Dx+(drhodx./rhoint).*funR).*uint + dudx.*funR + (Dy+(drhody./rhoint).*funR).*vint + dvdy.*funR ;

    for Inod=1:nod
        Jhh(:,:,Inod)=Jhh(:,:,Inod)+(W.*S(:,Inod)).*S ;
    end

end

nnzTot=nod*nod*Nele ; Iind=zeros(nnzTot,1) ; Jind=zeros(nnzTot,1) ; Xval=zeros(nnzTot,1) ; istak=0 ;
for Inod=1:nod
    for Knod=1:nod
        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod) ;
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod) ;
        Xval(istak+1:istak+Nele)=Jhh(:,Knod,Inod) ;
        istak=istak+Nele ;
    end
end
KJhh=sparse(Iind,Jind,Xval,Nnodes,Nnodes) ;

[~,dhdB,~,~,~,d2hdB2]=dGeometrydB(CtrlVar,F.s,F.S,F.B,F.b,F.rho,F.rhow) ;
[~,~,~,dhJhdot]=EvaluateJhdotAndDerivatives([],CtrlVar,MUA,F,[],Meas) ;      % BCs is not used in that function

Dh=spdiags(dhdB(:),0,Nnodes,Nnodes) ;
KhdotBB=Dh*KJhh*Dh+spdiags(dhJhdot(:).*d2hdB2(:),0,Nnodes,Nnodes) ;
KhdotBB=(KhdotBB+KhdotBB.')/2 ;

end
