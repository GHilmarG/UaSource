
function [UserVar,R,K]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,h0,ub0,vb0,a0,da0dh0,rho,h1,ub1,vb1,a1,da1dh1)

% Assembly



narginchk(15,15)
nargoutchk(2,3)

nOut=nargout;

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.hTheta;
dt=CtrlVar.dt;
CtrlVar.Tracer.SUPG.tau=CtrlVar.hSUPGtau;

h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(h1(MUA.connectivity,1),MUA.Nele,MUA.nod);

a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


Khh=zeros(MUA.Nele,MUA.nod,MUA.nod);
Rh=zeros(MUA.Nele,MUA.nod);


l=sqrt(2*MUA.EleAreas);
% vector over all elements for each  integration point

for Iint=1:MUA.nip  %Integration points
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    h0int=h0nod*fun;
    h1int=h1nod*fun;
    
    u0int=u0nod*fun; v0int=v0nod*fun;
    u1int=u1nod*fun; v1int=v1nod*fun;
    a0int=a0nod*fun;
    a1int=a1nod*fun;
    
    
    
    % derivatives at one integration point for all elements
    Deriv1=squeeze(Deriv(:,1,:)) ;
    Deriv2=squeeze(Deriv(:,2,:)) ;
    
    
    exx0=zeros(MUA.Nele,1);
    eyy0=zeros(MUA.Nele,1);
    
    exx1=zeros(MUA.Nele,1);
    eyy1=zeros(MUA.Nele,1);
    
    drhodx=zeros(MUA.Nele,1); drhody=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dhdx=dhdx+Deriv1(:,Inod).*h1nod(:,Inod);
        dhdy=dhdy+Deriv2(:,Inod).*h1nod(:,Inod);
        dh0dx=dh0dx+Deriv1(:,Inod).*h0nod(:,Inod);
        dh0dy=dh0dy+Deriv2(:,Inod).*h0nod(:,Inod);
        
        exx0=exx0+Deriv1(:,Inod).*u0nod(:,Inod);  % exx0
        eyy0=eyy0+Deriv2(:,Inod).*v0nod(:,Inod);
        
        
        drhodx=drhodx+Deriv1(:,Inod).*rhonod(:,Inod);
        drhody=drhody+Deriv2(:,Inod).*rhonod(:,Inod);
        
        exx1=exx1+Deriv1(:,Inod).*u1nod(:,Inod);
        eyy1=eyy1+Deriv2(:,Inod).*v1nod(:,Inod);
        
        
        
    end
    
    
    detJw=detJ*MUA.weights(Iint);
    
    
    speed0=sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.h.SUPG.tau) ;
    tauSUPGint=CtrlVar.SUPG.beta0*tau;
    
    
    
    qx1dx=rhoint.*exx.*h1int+rhoint.*u1int.*dh1dx+drhodx.*u1int.*hint;
    qy1dy=rhoint.*eyy.*h1int+rhoint.*v1int.*dh1dy+drhody.*v1int.*hint;
    qx0dx=rhoint.*exx0.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*h0int;
    qy0dy=rhoint.*eyy0.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*h0int;
    
    
    
    for Inod=1:MUA.nod
        
        
        SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
        
        
        if nOut>2
            for Jnod=1:MUA.nod
                
                Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                    +(rhoint.*fun(Jnod)...
                    -dt*theta*rhoint.*dadhint.*fun(Jnod)...
                    +dt*theta.*(rhoint.*exx.*fun(Jnod)+drhodx.*uint.*fun(Jnod)+rhoint.*uint.*Deriv(:,1,Jnod)+...
                    rhoint.*eyy.*fun(Jnod)+drhody.*vint.*fun(Jnod)+rhoint.*vint.*Deriv(:,2,Jnod)))...
                    .*SUPG.*detJw;
                
                
            end
        end
        
        % Note, I solve: LSH  \phi  = - RHS
        qterm=  dt*(theta*q1xdx+(1-theta)*q0xdx+theta*q1ydy+(1-theta)*q0ydy);
        dhdt=  rhoint.*(h0int-hint);
        accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int);
        
        rh= accterm - dhdt - qterm;
        
        % I solve K h = -R
        
        %
        Rh(:,Inod)=Rh(:,Inod)+rh.*SUPG.*detJw;
        
        
    end
end
%% assemble right-hand side

Rint=sparseUA(neq,1); Fext=sparseUA(neq,1);

for Inod=1:MUA.nod
    Rint=Rint+sparseUA(MUA.connectivity(:,Inod)+2*neqx,ones(MUA.Nele,1),Rh(:,Inod),neq,1);
end
%%

Rh=Rint-Fext;


if nargout>2
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

end