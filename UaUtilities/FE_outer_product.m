
function [UserVar,ab]=FE_outer_product(UserVar,CtrlVar,MUA,a,b)


%
% $$  a(x,y) b(x,y) $$
%
% $$ ( a_r \phi_r   \phi_p , \phi_q  \b_s \phi_s )   $$
%


narginchk(5,5)

ndim=2; dof=1; neq=dof*MUA.Nnodes;

a=a+zeros(MUA.Nnodes,1);
b=b+zeros(MUA.Nnodes,1);
anod=reshape(a(MUA.connectivity,1),MUA.Nele,MUA.nod);
bnod=reshape(b(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);


ElementMatrix=zeros(MUA.Nele,MUA.nod,MUA.nod);


% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    % values at integration point
    
    aint=anod*fun;
    bint=bnod*fun;
    
    detJw=detJ*weights(Iint);
    
    for Inod=1:MUA.nod
        
        
        for Jnod=1:MUA.nod

            ElementMatrix(:,Inod,Jnod)=ElementMatrix(:,Inod,Jnod)+...
                aint.*bint.*fun(Jnod).*fun(Inod).*detJw;
            
        end

    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=ElementMatrix(:,Inod,Jnod);
        istak=istak+MUA.Nele;
    end
end

ab=sparseUA(Iind,Jind,Xval,neq,neq);

end