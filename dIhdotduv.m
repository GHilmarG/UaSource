

function [UserVar,rh]=dIhdotduv(UserVar,CtrlVar,MUA,F,dhdtres,dhdtErr)

%
%  (hdot-hmeas) ( d(h du)/dx + d(h dv)/dv )
%

ndim=2; dof=2; neq=dof*MUA.Nnodes;
neqx=MUA.Nnodes ;


dhdtresnod=reshape(dhdtres(MUA.connectivity,1),MUA.Nele,MUA.nod);
dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),MUA.Nele,MUA.nod);

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);



% [points,weights]=sample('triangle',MUA.nip,ndim);

bx=zeros(MUA.Nele,MUA.nod);
by=zeros(MUA.Nele,MUA.nod);

% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration point
    
    dhdtresint=dhdtresnod*fun;
    dhdtErrint=dhdtErrnod*fun;
   
    hint=hnod*fun;
    
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    
    for Inod=1:MUA.nod
        
        termx=-dhdtresint.*(dhdx.*fun(Inod)+hint.*Deriv(:,1,Inod)).*detJw./dhdtErrint;
        termy=-dhdtresint.*(dhdy.*fun(Inod)+hint.*Deriv(:,2,Inod)).*detJw./dhdtErrint;
        
        bx(:,Inod)=bx(:,Inod)+termx;
        by(:,Inod)=by(:,Inod)+termy;
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),bx(:,Inod),neq,1);
    rh=rh+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),by(:,Inod),neq,1);
end

rh=full(rh) ; % This is a vector and I know it will not be particularly sparse


end