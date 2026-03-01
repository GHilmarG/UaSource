

function [dFhdotdu,dFhdotdv,dFhdotdhdot]=dFhdot_duvhdot(CtrlVar,MUA,F)




%%  Provides derivatives with respect to $u$, $v$, and $\dot{h}$ of $F_{\dot{h}}$
%
% $$F_{\dot{h}} = \dot{h} + ( \partial_x (u h ) + \partial (v h) ) -a = 0 $$
%
% The FE formulation is
%
% $$ \langle F_{\dot{h}} \, |\, \phi_k \rangle = 0 $$
%
% for $k=1\ldots n$
%
% And the $u$ derivataive is
%
% $$ \langle \partial_u F_{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle $$
%
% where
%
% $$ \partial_u \dot{h} = \partial_ u  (  a - ( \partial_x (u h ) + \partial (v h) )) = - \partial_ u  \,  \partial_x (u h ) $$
%
% or
%
% $$ \partial_u \dot{h} = - \partial_x  \,  \partial_u (u h ) = -\partial_x (\delta u \, h) = - h \, \partial_x \delta u - \delta u \, \partial_x h $$
%
%
% Therefore:
%
% $$
% \langle \partial_u F_{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle - h \, \partial_x \phi_i - \phi_i \, \partial_x h | \phi_k \rangle
% $$
%
% $$
% \langle \partial_v F_{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle - h \, \partial_y \phi_i - \phi_i \, \partial_y h | \phi_k \rangle
% $$
%
% $$
% \langle \partial_{\dot{h}}  F_{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle \phi_i | \phi_ k \rangle
% $$
%
% see also: dhdtExplicit.m
%
%%

ndim=2;
nNodes=MUA.Nnodes ;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

dFhdotdu=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFhdotdv=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFhdotdhdot=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    h=hnod*fun;

    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dudx=dudx+Deriv(:,1,Inod).*ubnod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vbnod(:,Inod);

        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);


    end

    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            dFhdotdu(:,Inod,Jnod)=dFhdotdu(:,Inod,Jnod) - (h.* Deriv(:,1,Inod) + fun(Jnod).*dhdx )   .*detJw;
            dFhdotdv(:,Inod,Jnod)=dFhdotdv(:,Inod,Jnod) - (h.* Deriv(:,2,Inod) + fun(Jnod).*dhdy )   .*detJw;
            dFhdotdhdot(:,Inod,Jnod)=dFhdotdhdot(:,Inod,Jnod) + fun(Inod).*fun(Jnod) .*detJw;

        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');

dFhdotduXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
dFhdotdvXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
dFhdotdhdotXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);

        dFhdotduXval(istak+1:istak+MUA.Nele)=dFhdotdu(:,Inod,Jnod);
        dFhdotdvXval(istak+1:istak+MUA.Nele)=dFhdotdv(:,Inod,Jnod);
        dFhdotdhdotXval(istak+1:istak+MUA.Nele)=dFhdotdhdot(:,Inod,Jnod);


    end
end


dFhdotdu=sparseUA(Iind,Jind,dFhdotduXval,2*nNodes,nNodes);
dFhdotdv=sparseUA(Iind,Jind,dFhdotdvXval,2*nNodes,nNodes);
dFhdotdhdot=sparseUA(Iind,Jind,dFhdotdhdotXval,2*nNodes,nNodes);

end








