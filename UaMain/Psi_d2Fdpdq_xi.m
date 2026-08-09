


function [Psi_d2FdCC]=Psi_d2Fdpdq_xi(CtrlVar,MUA,F,uAdjoint,vAdjoint,KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC)

narginchk(14,14)

%%
%
% Calculates
%
% $$  \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj} $$
%
% and
%
% $$ \Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki} $$
%
%
% which is a part of the Hessian.
%
% This term also involves mixed derivatives of the forward model, and is not zero
%
%
% For
%
% $$ H^{\#8}_{ij}=\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki} $$
%
% and
%
% $$ H^{\#6}_{ij} = \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj} $$
%
% we have
%
% $$ [H^{\#8}]^T=H^{\#8}_{ji}=\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_i} \xi_{kj} =\Psi_n \frac{\partial^2 F_n}{\partial p_i \,\partial q_k} \xi_{kj} = H^{\#6}_{ij}= [H^{\#6}] $$
%
% So these two terms are just the transpose of each other, and both are calculated here.
%
% Also, here $q$ is both $u$ and $v$
%
% $$ H^{\#8}_{ij}=
% \Psi^x_n \frac{\partial^2 F^x_n}{ \partial C_i \, \partial u_k} \frac{\partial u_k}{\partial C_j} + \Psi^x_n \frac{\partial^2 F^x_n}{ \partial C_i \, \partial v_k} \frac{\partial v_k}{\partial C_j}
% +\Psi^y_n \frac{\partial^2 F^y_n}{ \partial C_i \, \partial u_k} \frac{\partial u_k}{\partial C_j} + \Psi^y_n \frac{\partial^2 F^y_n}{ \partial C_i \, \partial v_k} \frac{\partial v_k}{\partial C_j}
% $$
%
%
%
% $$ \langle \Psi(x) \mid \delta^2_{uC} F(x) \rangle = \langle \Psi(x) \mid \partial^2_{uC} F(x) \, \delta u \, \delta C \rangle $$
%
% Or in a discrete form:
%
% $$\Psi_n \frac{\partial^2 F_n}{\partial u_i \, \partial C_j} $$
%
% For example using Weertman sliding law:
%
% $$ F^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi_i \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, u | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% $$ F^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi_i \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, v | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
% where
%
% $$\beta^2= (C+C_0)^{-1/m} \, (u^2+v^2+v_0^2)^{(1/m-1)/2} $$
%
% the relevant part is
%
% $$F_x= \mathcal{G} \beta^2 \, v_x $$
%
% $$F_y= \mathcal{G} \beta^2 \, v_y $$
%
%
% $$ \beta^2 = (C+C_0)^{-1/m} \, (v_x^2+v_y^2 + v_0)^{(1-m)/2m} $$
%
% $$F_x= \mathcal{G} \; (C+C_0)^{-1/m} \;  (v_x^2+v_y^2+v_0^2)^{(1-m)/2m}  \, v_x $$
%
% we have
%
% $$\delta_C F_x = -\frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; (v_x^2+v_y^2 + v_0^2)^{(1-m)/2m} \, v_x \; \delta C$$
%
% and
%
% $$\delta^2_{v_x \,C } F_x = -\frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; \left  (\frac{1-m}{m} (v_x^2+v_y^2 + v_0^2)^{(1-3m)/2m} \, v_x \, v_x +  (v_x^2+v_y^2 + v_0)^{(1-m)/2m} \right )    \; \delta C \, \delta u$$
%
%
%
% $$ \langle \Psi_x | \delta^2_{v_x C} F_x \rangle =
%  -\int  \Psi_x \,  \frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; \left  (\frac{1-m}{m} (v_x^2+v_y^2 + v_0)^{(1-3m)/2m} \, v_x \, v_x +  (v_x^2+v_y^2 + v_0^2)^{(1-m)/2m} \right )    \; \phi_i \, \phi_j \; dx \, dy
% $$
%
% which is a $n \times n$ matrix, which then needs to be multiplied with the $n \times n$ sensitivity matrix $\partial
% v_x/\partial C$
%
% There are two components to $F=(F^x,F^y)$ and two components to $q=(v_x,v_y)$, so we get
%
% $$ \langle \Psi_x | \delta^2_{v_x C} F_x \rangle \, \delta_c v_x + \langle  \Psi_x | \delta^2_{v_y C} F_x   \rangle \, \delta_C v_y   
% +  \langle \Psi_y | \delta^2_{v_x C} F_y \rangle \, \delta_C v_x + \langle  \Psi_y | \delta^2_{v_y C} F_y   \rangle \, \delta_C v_y $$
%
%
%%

if CtrlVar.SlidingLaw~="Weertman"

    error("Psi_d2Fdpdq_xi:NotImplemented","only implemented for Weertman sliding law.")

end

if ~contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)

    error("Psi_d2Fdpdq_xi:NotImplemented","only implemented for log(C) .")

end


ndim=2;
nNodes=MUA.Nnodes ;

C0=CtrlVar.Cmin;
v0=CtrlVar.SpeedZero;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod

Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hfnod=F.rhow*(Snod-Bnod)./rhonod;

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);



Hu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Hv=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    % Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    h=hnod*fun;
    hf=hfnod*fun;

    u=unod*fun;
    v=vnod*fun;



    C=Cnod*fun;
    C(C<C0)=C0;

    m=mnod*fun;

    lx=uAdjointnod*fun;
    ly=vAdjointnod*fun;

    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

    detJw=detJ*MUA.weights(Iint);


    
    % dcFx=G.*(-1./m).* (C+C0).^(-1./m-1) .*  (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)).*u ;
    % dcFy=G.*(-1./m).* (C+C0).^(-1./m-1) .*  (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)).*v ;


    ducFx=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)-1) .*u.*u + (u.^2+v.^2 + v0^2).^((1-m)./(2.*m))  ) ;
    dvcFx=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)-1) .*v.*u  ) ;

    ducFy=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)-1) .*u.*v  ) ;
    dvcFy=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)-1) .*v.*v + (u.^2+v.^2 + v0^2).^((1-m)./(2.*m))  ) ;



    l_d2FdudC=(lx.*ducFx+ly.*ducFy).*C.*log(10);  % There is only one derivative with respect to C, so chain rule is only used once
    l_d2FdvdC=(lx.*dvcFx+ly.*dvcFy).*C.*log(10);


    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            Hu(:,Inod,Jnod)=Hu(:,Inod,Jnod) + l_d2FdudC .*fun(Inod) .*fun(Jnod).*detJw;
            Hv(:,Inod,Jnod)=Hv(:,Inod,Jnod) + l_d2FdvdC .*fun(Inod) .*fun(Jnod).*detJw;

        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
HuVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
HvVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        HuVal(istak+1:istak+MUA.Nele)=Hu(:,Inod,Jnod);
        HvVal(istak+1:istak+MUA.Nele)=Hv(:,Inod,Jnod);


        istak=istak+MUA.Nele;

    end
end

Psid2FdCdu=sparseUA(Iind,Jind,HuVal,nNodes,nNodes) ;  % this will be full matrix,
Psid2FdCdv=sparseUA(Iind,Jind,HuVal,nNodes,nNodes) ;  % this will be full matrix,


Psi_d2FduC_dudC=Psid2FdCdu*KdudC;  % this will be full matrix,
Psi_d2FdvC_dvdC=Psid2FdCdv*KdvdC;

Psi_d2FduC_dudC=full(Psi_d2FduC_dudC);
Psi_d2FdvC_dvdC=full(Psi_d2FdvC_dvdC);

Psi_d2FduC_dudC=Psi_d2FduC_dudC+Psi_d2FduC_dudC';
Psi_d2FdvC_dvdC=Psi_d2FdvC_dvdC+Psi_d2FdvC_dvdC';

Psi_d2FdCC=Psi_d2FduC_dudC+Psi_d2FdvC_dvdC;


end







