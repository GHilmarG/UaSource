

function [UserVar,kv,rh]=TracerConservationEquationAssembly(UserVar,CtrlVar,MUA,dt,h0,u0,v0,a0,u1,v1,a1,kappa)

% Note: Assembly for the linear tracer equation
%  dc/dt + d (u c)/dx + d (v c)/dy - div (kappa grad c) = a

%%
%
% $$ \partial c/\partial t + \partial  ( u c)/ \partial x + \partial  ( v c)/\partial y - \nabla \cdot (\kappa \nabla c ) = a$$
%
%
%    
%%

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.theta;
%tauSUPG=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0,v0,dt,MUA);


h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);
v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);
a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);

if numel(kappa)==1
    kappa=kappa+zeros(MUA.Nnodes,1);
end

kappanod=reshape(kappa(MUA.connectivity,1),MUA.Nele,MUA.nod);
%tauSUPGnod=reshape(tauSUPG(MUA.connectivity,1),MUA.Nele,MUA.nod);


% [points,weights]=sample('triangle',MUA.nip,ndim);


d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
b1=zeros(MUA.Nele,MUA.nod);

l=sqrt(2*MUA.EleAreas);

% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    
 

    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    
    
    
    
    
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration point
    
    h0int=h0nod*fun;
    u0int=u0nod*fun;
    v0int=v0nod*fun;
    a0int=a0nod*fun;
    
    u1int=u1nod*fun;
    v1int=v1nod*fun;
    a1int=a1nod*fun;
    

    
    
    kappaint=kappanod*fun;
    %tauSUPGint=tauSUPGnod*fun;
    
    du1dx=zeros(MUA.Nele,1); du0dx=zeros(MUA.Nele,1); dh0dx=zeros(MUA.Nele,1);
    dv1dy=zeros(MUA.Nele,1); dv0dy=zeros(MUA.Nele,1); dh0dy=zeros(MUA.Nele,1);
    
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        du1dx=du1dx+Deriv(:,1,Inod).*u1nod(:,Inod);
        du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
        
        dv1dy=dv1dy+Deriv(:,2,Inod).*v1nod(:,Inod);
        dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
        
        dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
        dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    
    % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
    %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
    
    
    % SUPG parameter calculation taken inside int-loop on 20 Sept, 2021
    
    speed0=sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.Tracer.SUPG.tau) ; 
    tauSUPGint=CtrlVar.SUPG.beta0*tau;
    
   
    
    
    for Inod=1:MUA.nod
        
        SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
        
        SUPGdetJw=SUPG.*detJw;
        
        for Jnod=1:MUA.nod
            
            h1term=fun(Jnod).*SUPGdetJw;
            
            
                  
            hdxu1=dt*theta*du1dx.*fun(Jnod).*SUPGdetJw;
            udxh1=dt*theta*u1int.*Deriv(:,1,Jnod).*SUPGdetJw;
            
            hdyv1=dt*theta*dv1dy.*fun(Jnod).*SUPGdetJw;
            vdyh1=dt*theta*v1int.*Deriv(:,2,Jnod).*SUPGdetJw;
            
            kdxh1=dt*theta.*kappaint.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
            kdyh1=dt*theta.*kappaint.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;
            
            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+h1term+hdxu1+udxh1+hdyv1+vdyh1+kdxh1+kdyh1;
            
        end
        
        h0term=h0int.*SUPGdetJw;  % this is the h term, ie not \Delta h term (because the system is linear in h no need to write it in incremental form)
        a0term=dt*(1-theta)*a0int.*SUPGdetJw;
        a1term=dt*theta*a1int.*SUPGdetJw;
        
        hdxu0=-dt*(1-theta)*du0dx.*h0int.*SUPGdetJw;
        udxh0=-dt*(1-theta)*dh0dx.*u0int.*SUPGdetJw;
        
        hdyv0=-dt*(1-theta)*dv0dy.*h0int.*SUPGdetJw;
        vdyh0=-dt*(1-theta)*dh0dy.*v0int.*SUPGdetJw;
        
        kappadhdx=-dt*(1-theta)*kappaint.*dh0dx.*Deriv(:,1,Inod).*detJw;
        kappadhdy=-dt*(1-theta)*kappaint.*dh0dy.*Deriv(:,2,Inod).*detJw;
        
        b1(:,Inod)=b1(:,Inod)+h0term+a0term+a1term+hdxu0+udxh0+hdyv0+vdyh0+kappadhdx+kappadhdy;
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b1(:,Inod),neq,1);
end

% tic
% kv=sparse(neq,neq);
% for Inod=1:nod
%     for Jnod=1:nod
%         kv=kv+sparseUA(MUA.connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
%     end
% end
% toc


Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
        istak=istak+MUA.Nele;
    end
end

kv=sparseUA(Iind,Jind,Xval,neq,neq);






end