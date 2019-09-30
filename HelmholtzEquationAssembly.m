

function [UserVar,Matrix,rh]=HelmholtzEquationAssembly(UserVar,CtrlVar,MUA,a,b,c,d)

%
%  $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y)$$
%

narginchk(7,7)


ndim=2; dof=1; neq=dof*MUA.Nnodes;

a=a+zeros(MUA.Nnodes,1);
b=b+zeros(MUA.Nnodes,1);
c=c+zeros(MUA.Nnodes,1);
d=d+zeros(MUA.Nnodes,1);

anod=reshape(a(MUA.connectivity,1),MUA.Nele,MUA.nod);
bnod=reshape(b(MUA.connectivity,1),MUA.Nele,MUA.nod);
cnod=reshape(c(MUA.connectivity,1),MUA.Nele,MUA.nod);
dnod=reshape(d(MUA.connectivity,1),MUA.Nele,MUA.nod);


[points,weights]=sample('triangle',MUA.nip,ndim);


ElementMatrix=zeros(MUA.Nele,MUA.nod,MUA.nod);
ElementRHS=zeros(MUA.Nele,MUA.nod);

% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration point
    
    aint=anod*fun;
    bint=bnod*fun;
    cint=cnod*fun;
    
    dddx=zeros(MUA.Nele,1); 
    dddy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dddx=dddx+Deriv(:,1,Inod).*dnod(:,Inod);
        dddy=dddy+Deriv(:,2,Inod).*dnod(:,Inod);
                
    end
    
    
    detJw=detJ*weights(Iint);
    
    % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
    %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
    
    
    
    for Inod=1:MUA.nod
        
        
        for Jnod=1:MUA.nod
            
            aphiphi=aint.*fun(Jnod).*fun(Inod).*detJw;
            
            bdphidxdphidx=bint.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
            bdphidydphidy=bint.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;
            
            ElementMatrix(:,Inod,Jnod)=ElementMatrix(:,Inod,Jnod)+aphiphi+bdphidxdphidx+bdphidydphidy ;
            
        end
        
        ElementRHS(:,Inod)=ElementRHS(:,Inod)+...
            cint.*fun(Inod).*detJw + ...
            dddx.*Deriv(:,1,Inod).*detJw + ...
            dddy.*Deriv(:,2,Inod).*detJw ;
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),ElementRHS(:,Inod),neq,1);
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

Matrix=sparseUA(Iind,Jind,Xval,neq,neq);


end