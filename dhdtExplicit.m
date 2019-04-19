
function [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F)

%
%
%  dh/dt=a -  ( dqx/dx + dqy/dy)
%
%  (hdot-hmeas) ( d(h du)/dx + d(h dv)/dv )
%

ndim=2; dof=1; neq=dof*MUA.Nnodes;



anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);

%qx=F.h.*F.ub;
%qy=F.h.*F.vb;
%qxnod=reshape(qx(MUA.connectivity,1),MUA.Nele,MUA.nod);
%qynod=reshape(qy(MUA.connectivity,1),MUA.Nele,MUA.nod);

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);


[points,weights]=sample('triangle',MUA.nip,ndim);

b=zeros(MUA.Nele,MUA.nod);


% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    aint=anod*fun;
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        %dqxdx=dqxdx+Deriv(:,1,Inod).*qxnod(:,Inod);
        %dqydy=dqydy+Deriv(:,2,Inod).*qynod(:,Inod);
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
    end
    
    detJw=detJ*weights(Iint);
    
    for Inod=1:MUA.nod
        
        tx=(dhdx.*uint+hint.*dudx);
        ty=(dhdy.*vint+hint.*dvdy);
        
        term=(aint-tx-ty).*fun(Inod).*detJw;
        b(:,Inod)=b(:,Inod)+term;
        
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b(:,Inod),neq,1);
end

dhdt=MUA.M\rh;
dhdt=full(dhdt);

end