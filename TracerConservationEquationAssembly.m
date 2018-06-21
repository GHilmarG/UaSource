

function [UserVar,kv,rh]=TracerConservationEquationAssembly(UserVar,CtrlVar,MUA,dt,h0,u0,v0,a0,u1,v1,a1,kappa)

coordinates=MUA.coordinates;
connectivity=MUA.connectivity; 
nip=MUA.nip; 

Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
ndim=2; dof=1; neq=dof*Nnodes;

theta=CtrlVar.theta;

h0nod=reshape(h0(connectivity,1),Nele,nod);
u0nod=reshape(u0(connectivity,1),Nele,nod);   % Nele x nod
u1nod=reshape(u1(connectivity,1),Nele,nod);
v0nod=reshape(v0(connectivity,1),Nele,nod);   % Nele x nod
v1nod=reshape(v1(connectivity,1),Nele,nod);
a0nod=reshape(a0(connectivity,1),Nele,nod);
a1nod=reshape(a1(connectivity,1),Nele,nod);

if numel(kappa)==1
    kappa=kappa+zeros(MUA.Nnodes,1);
end

kappanod=reshape(kappa(connectivity,1),Nele,nod);


[points,weights]=sample('triangle',nip,ndim);


d1d1=zeros(Nele,nod,nod);
b1=zeros(Nele,nod);


% SUPG specific

l=2*sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));



% vector over all elements for each integartion point
for Iint=1:nip
    
    
    
    fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
   
    
    
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
    
    du1dx=zeros(Nele,1); du0dx=zeros(Nele,1); dh0dx=zeros(Nele,1);
    dv1dy=zeros(Nele,1); dv0dy=zeros(Nele,1); dh0dy=zeros(Nele,1);
    
    % derivatives at one integration point for all elements
    for Inod=1:nod
        
        du1dx=du1dx+Deriv(:,1,Inod).*u1nod(:,Inod);
        du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
        
        dv1dy=dv1dy+Deriv(:,2,Inod).*v1nod(:,Inod);
        dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
        
        dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
        dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
        
    end
    
    detJw=detJ*weights(Iint);
    
    % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
    %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
    
    % SUPG weight
    
    speed0=sqrt(u0int.*u0int+v0int.*v0int+eps^2);
    speed1=sqrt(u1int.*u1int+v1int.*v1int+eps^2);
    
    ECN=speed0.*dt./l+eps; % This is the `Element Courant Number' ECN
    
    K=coth(ECN)-1./ECN;
    
    
    
    tau0=CtrlVar.SUPG.beta0*K.*l./speed0 ; % sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
    tau1=CtrlVar.SUPG.beta1*K.*l./speed1 ; % sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);
    
    
    % SUPG introduces diffusion with a diffusion coefficient of
    % u1int*SUPG = u1int* tau*u1int
    %            = u1int^2*length/speed
    %             with units: length^2/time
    %
    % kappaSUPG=beta*K*u^2*l/speed
    %
    % kappaSUPGu=CtrlVar.SUPG.beta0*K.*u0int.^2.*l./speed0 ;
    %
    % However, because the SUPG weighting term is applied to all terms
    % including the time derivative term (ie. consistent formulation) 
    % the net results is not just that of a single added diffusion term. 
    %    
    
    for Inod=1:nod
        
        
        SUPG=fun(Inod)+theta.*tau1.*(u1int.*Deriv(:,1,Inod)+v1int.*Deriv(:,2,Inod))...
            +(1-theta).*tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));

        
        for Jnod=1:nod
            
            h1term=fun(Jnod).*SUPG.*detJw;
            
            hdxu1=dt*theta*du1dx.*fun(Jnod).*SUPG.*detJw;
            udxh1=dt*theta*u1int.*Deriv(:,1,Jnod).*SUPG.*detJw;
            
            hdyv1=dt*theta*dv1dy.*fun(Jnod).*SUPG.*detJw;
            vdyh1=dt*theta*v1int.*Deriv(:,2,Jnod).*SUPG.*detJw;
            
            kdxh1=dt*theta.*kappaint.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
            kdyh1=dt*theta.*kappaint.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;
            
            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+h1term+hdxu1+udxh1+hdyv1+vdyh1+kdxh1+kdyh1;

        end
        
        h0term=h0int.*SUPG.*detJw;  % this is the h term, ie not \Delta h term (because the system is linear in h no need to write it in incremental form)
        a0term=dt*(1-theta)*a0int.*SUPG.*detJw;
        a1term=dt*theta*a1int.*SUPG.*detJw;
        
        hdxu0=-dt*(1-theta)*du0dx.*h0int.*SUPG.*detJw;
        udxh0=-dt*(1-theta)*dh0dx.*u0int.*SUPG.*detJw;
        
        hdyv0=-dt*(1-theta)*dv0dy.*h0int.*SUPG.*detJw;
        vdyh0=-dt*(1-theta)*dh0dy.*v0int.*SUPG.*detJw;
        
        kappadhdx=-dt*(1-theta)*kappaint.*dh0dx.*Deriv(:,1,Inod).*detJw;
        kappadhdy=-dt*(1-theta)*kappaint.*dh0dy.*Deriv(:,2,Inod).*detJw;
        
        b1(:,Inod)=b1(:,Inod)+h0term+a0term+a1term+hdxu0+udxh0+hdyv0+vdyh0+kappadhdx+kappadhdy;
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:nod
    rh=rh+sparseUA(connectivity(:,Inod),ones(Nele,1),b1(:,Inod),neq,1);
end

% tic
% kv=sparse(neq,neq);
% for Inod=1:nod
%     for Jnod=1:nod
%         kv=kv+sparseUA(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
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