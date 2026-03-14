

function [KdFhdotduvhdot]=dFhdot_duvhdot(CtrlVar,MUA,F)

narginchk(3,3)
nargoutchk(1,1)


%%  Provides derivatives of $F^{\dot{h}}$ with respect to $u$, $v$, and $\dot{h}$.
%
% $$
% \left (\begin{array}{ccc}
% \frac{\partial F^{\dot{h}}}{\partial u}  & \frac{\partial F^{\dot{h}}}{\partial v}  & \frac{\partial F^{\dot{h}}}{\partial h}
% \end{array}\right )
% $$
%
% This is a $n \times 3n$ matrix
%
%
% $$F^{\dot{h}} = \dot{h} + ( \partial_x (u h ) + \partial (v h) ) -a = 0 $$
%
% The FE formulation is
%
% $$ \langle F^{\dot{h}} \, |\, \phi_k \rangle = 0 $$
%
% for $k=1\ldots n$
%
% And the $u$ derivative is
%
% $$ \langle \partial_u F{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle $$
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
% \langle \partial_u F^{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle -h \, \partial_x \phi_i + \phi_i \, \partial_x h | \phi_k \rangle
% $$
%
% $$
% \langle \partial_v F^{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle  h \, \partial_y \phi_i + \phi_i \, \partial_y h | \phi_k \rangle
% $$
%
% $$
% \langle \partial^{\dot{h}}  F^{\dot{h}} \, \phi_i \, | \, \phi_k  \rangle = \langle \phi_i | \phi_ k \rangle
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

            dFhdotdu(:,Inod,Jnod)=dFhdotdu(:,Inod,Jnod) + (h.* Deriv(:,1,Jnod) + fun(Jnod).*dhdx).*fun(Inod).*detJw;
            dFhdotdv(:,Inod,Jnod)=dFhdotdv(:,Inod,Jnod) + (h.* Deriv(:,2,Jnod) + fun(Jnod).*dhdy).*fun(Inod) .*detJw;
            dFhdotdhdot(:,Inod,Jnod)=dFhdotdhdot(:,Inod,Jnod) + fun(Jnod).*fun(Inod) .*detJw;  % I don't really need to do this as this is the Mass Matrix

        end
    end
end


%% This is if I want three matrices as an output
% Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
% Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
% 
% dFhdotduXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
% dFhdotdvXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
% dFhdotdhdotXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
% 
% istak=0;
% 
% for Inod=1:MUA.nod
%     %istak=0;
% 
%     for Jnod=1:MUA.nod
% 
%         Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
%         Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
% 
%         dFhdotduXval(istak+1:istak+MUA.Nele)=dFhdotdu(:,Inod,Jnod);
%         dFhdotdvXval(istak+1:istak+MUA.Nele)=dFhdotdv(:,Inod,Jnod);
%         dFhdotdhdotXval(istak+1:istak+MUA.Nele)=dFhdotdhdot(:,Inod,Jnod);
% 
%         istak=istak+MUA.Nele;
% 
%     end
% end
% 
% 
% KdFhdotdu=sparseUA(Iind,Jind,dFhdotduXval,nNodes,nNodes);
% KdFhdotdv=sparseUA(Iind,Jind,dFhdotdvXval,nNodes,nNodes);
% KdFhdotdhdot=sparseUA(Iind,Jind,dFhdotdhdotXval,nNodes,nNodes);
% 



%% This creates one $n \times 3n$ matrix
Iind=zeros(MUA.nod*MUA.nod*3*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*3*MUA.Nele,1,'uint32');
dFhduvhXval=zeros(MUA.nod*MUA.nod*3*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        dFhduvhXval(istak+1:istak+MUA.Nele)=dFhdotdu(:,Inod,Jnod);
        
        
        istak=istak+MUA.Nele;
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+nNodes; 
        dFhduvhXval(istak+1:istak+MUA.Nele)=dFhdotdv(:,Inod,Jnod);

     
        istak=istak+MUA.Nele;
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*nNodes; 
        dFhduvhXval(istak+1:istak+MUA.Nele)=dFhdotdhdot(:,Inod,Jnod);

        istak=istak+MUA.Nele;

    end
end

KdFhdotduvhdot=sparseUA(Iind,Jind,dFhduvhXval,nNodes,3*nNodes);



end








