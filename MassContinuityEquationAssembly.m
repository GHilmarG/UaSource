
function [UserVar,R,K]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,h0,ub0,vb0,a0,da0dh0,rho0,h1,ub1,vb1,a1,da1dh1,rho1)

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

q0x=rho0.*h0.*ub0 ;  q0y=rho0.*h0.*vb0; 
q1x=rho1.*h1.*ub1 ;  q1y=rho1.*h1.*vb1; 


q0xnod=reshape(q0x(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
q0ynod=reshape(q0y(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
q1xnod=reshape(q1x(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
q1ynod=reshape(q1y(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod



a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);


d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
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

    dq0xdx=zeros(MUA.Nele,1); dq0ydy=zeros(MUA.Nele,1);
    dq1xdx=zeros(MUA.Nele,1); dq1ydy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dq0xdx=dq0xdx+Deriv(:,1,Inod).*q0xnod(:,Inod);
        dq0ydy=dq0ydy+Deriv(:,2,Inod).*q0ynod(:,Inod);
        dq1xdx=dq1xdx+Deriv(:,1,Inod).*q1xnod(:,Inod);
        dq1ydy=dq1ydy+Deriv(:,2,Inod).*q1ynod(:,Inod);
      
    end
    
        
    detJw=detJ*MUA.weights(Iint);
    
        
    speed0=sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.h.SUPG.tau) ;
    tauSUPGint=CtrlVar.SUPG.beta0*tau;
    
    nod=MUA.nod ;
    
    for Inod=1:MUA.nod
        
        
        SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
        
        
        if nOut>2
            for Jnod=1:nod
                
                
                TL=fun(Jnod);

      
                LL=dt*theta*(...
                    u1int.*Deriv(:,1,Jnod) + v1int.*Deriv(:,2,Jnod)...
                    -a1int.*(n1xDeriv1Jnod + n1yDeriv2Jnod)...
                    );
                
                % Pertubation term (diffusion)
       
                
                AddUp=((isT*TL+isL*LL).*(fun(Inod)+isPG.*SUPG)+isP*Plhs).*detJw ;

         
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+AddUp;
                
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