function [UserVar,R,K]=NEquationAssembly(UserVar,CtrlVar,MUA,F0,F1,K,Sw,Phi) 



%% Effective water pressure equation
%
% 
% $$\frac{S_w}{\rho g}  \partial_t N +  \mathbf{V} \cdot \nabla N  - \nabla \cdot (\kappa h_w \nabla N) = m + \nabla \cdot ( k h_w ((\rho_w-\rho) g \nabla b - \rho g \nabla s)) $$ 
%
% $$\frac{S_w}{\rho g}  \partial_t N +  \mathbf{V} \cdot \nabla N  - \nabla \cdot (\kappa h_w \nabla N) = m + \nabla \cdot (k h_w \Phi )  $$ 
%
% 
%%


ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.theta;
dt=CtrlVar.dt ; 
%tauSUPG=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0,v0,dt,MUA);

N0nod=reshape(F0.N(MUA.connectivity,1),MUA.Nele,MUA.nod);

s0nod=reshape(F0.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
s1nod=reshape(F1.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

b0nod=reshape(F0.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
b1nod=reshape(F1.b(MUA.connectivity,1),MUA.Nele,MUA.nod);

% Vx0nod=reshape(Vx0(MUA.connectivity,1),MUA.Nele,MUA.nod);   
% Vy0nod=reshape(Vy0(MUA.connectivity,1),MUA.Nele,MUA.nod);   
% Vx1nod=reshape(Vx1(MUA.connectivity,1),MUA.Nele,MUA.nod);   
% Vy1nod=reshape(Vy1(MUA.connectivity,1),MUA.Nele,MUA.nod);   
% 

aw0nod=reshape(F0.aw(MUA.connectivity,1),MUA.Nele,MUA.nod);
aw1nod=reshape(F1.aw(MUA.connectivity,1),MUA.Nele,MUA.nod);


Knod=reshape(K(MUA.connectivity,1),MUA.Nele,MUA.nod);


if isempty(Phi)
    Phi=F.g.* ((F.rhow-F.rho).*F.b + F.rho.* F.s);

end

Sw=Sw./(F0.rho.*F0.g) ; 
Swnod=reshape(Sw(MUA.connectivity,1),MUA.Nele,MUA.nod);
Phinod=reshape(Phi(MUA.connectivity,1),MUA.Nele,MUA.nod);





d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
b1=zeros(MUA.Nele,MUA.nod);

l=sqrt(2*MUA.EleAreas);

% vector over all elements for each integration point
for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);


    % Deriv : Nele x dof x nod
    %  detJ : Nele

    % values at integration point
    
    N0int=N0nod*fun;
    % 
    % Vx0int=Vx0nod*fun; Vy0int=Vy0nod*fun;
    % Vx1int=Vx1nod*fun; Vy1int=Vy1nod*fun;
    % 
    aw0int=aw0nod*fun;
    aw1int=aw1nod*fun;
    
    K=Knod*fun;
    Sw=Swnod*fun; 
    
    
    
    ds0dx=zeros(MUA.Nele,1); ds0dy=zeros(MUA.Nele,1); 
    ds1dx=zeros(MUA.Nele,1); ds1dy=zeros(MUA.Nele,1); 
    
    db0dx=zeros(MUA.Nele,1); db0dy=zeros(MUA.Nele,1); 
    db1dx=zeros(MUA.Nele,1); db1dy=zeros(MUA.Nele,1); 
    
    dN0dx=zeros(MUA.Nele,1); dN0dy=zeros(MUA.Nele,1); 
    dPhidx=zeros(MUA.Nele,1); dPhidy=zeros(MUA.Nele,1); 
    
    
 
    
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        ds0dx=ds0dx+Deriv(:,1,Inod).*s0nod(:,Inod);
        ds0dy=ds0dy+Deriv(:,2,Inod).*s0nod(:,Inod);

        db0dx=db0dx+Deriv(:,1,Inod).*b0nod(:,Inod);
        db0dy=db0dy+Deriv(:,2,Inod).*b0nod(:,Inod);

        ds1dx=ds1dx+Deriv(:,1,Inod).*s1nod(:,Inod);
        ds1dy=ds1dy+Deriv(:,2,Inod).*s1nod(:,Inod);

        db1dx=db1dx+Deriv(:,1,Inod).*b1nod(:,Inod);
        db1dy=db1dy+Deriv(:,2,Inod).*b1nod(:,Inod);
        
        dN0dx=dN0dx+Deriv(:,1,Inod).*N0nod(:,Inod);
        dN0dy=dN0dy+Deriv(:,2,Inod).*N0nod(:,Inod);

        dPhidx=dPhidx+Deriv(:,1,Inod).*Phinod(:,Inod);
        dPhidy=dPhidy+Deriv(:,2,Inod).*Phinod(:,Inod);
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    
    % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
    %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
    
    
    % SUPG parameter calculation taken inside int-loop on 20 Sept, 2021
    
    % speed0=sqrt(Vx0int.*Vx0int+Vy0int.*Vy0int+CtrlVar.SpeedZero^2);
    speed0=0 ; 

    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.Tracer.SUPG.tau) ; 
    tauSUPGint=CtrlVar.SUPG.beta0*tau;
    

    for Inod=1:MUA.nod

        % SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(Vx0int.*Deriv(:,1,Inod)+Vy0int.*Deriv(:,2,Inod));
        SUPG=fun(Inod) ;   % not using this for the time being

        SUPGdetJw=SUPG.*detJw;

        for Jnod=1:MUA.nod

            N1term=Sw.*fun(Jnod).*fun(Inod).*detJw ; % 

            %
            % Vx1dN1dx=dt*theta*Vx1int.*Deriv(:,1,Jnod).*SUPGdetJw;
            % Vy1dN1dy=dt*theta*Vy1int.*Deriv(:,2,Jnod).*SUPGdetJw;
            Vx1dN1dx=0 ;
            Vy1dN1dy=0 ;

            KdN1dx=-dt*theta.*K.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
            KdN1dy=-dt*theta.*K.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;

            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+N1term+Vx1dN1dx+Vy1dN1dy+KdN1dx+KdN1dy;

        end

        N0term=Sw.*N0int.*fun(Inod).*detJw;  % this is the h term, ie not \Delta h term (because the system is linear in h no need to write it in incremental form)

        % Vx0dN0dx=-dt*(1-theta)*Vx0int.*dN0dx.*SUPGdetJw;
        % Vy0dN0dy=-dt*(1-theta)*Vy0int.*dN0dy.*SUPGdetJw;
        Vx0dN0dx=0 ;
        Vy0dN0dy=0 ;

        NTerm=dt*(1-theta)*K.*(dN0dx.*Deriv(:,1,Inod)+dN0dy.*Deriv(:,2,Inod)) .*detJw;

        % total contribution, to be changed later
        % PhiTerm =-dt*Kint.*g.*  (((rhow-rho).*db0dx-rho.*ds0dx) .*Deriv(:,1,Inod) + ((rhow-rho).*db0dy-rho.*ds0dy) .*Deriv(:,2,Inod) ).*detJw;

        PhiTerm=-dt*K.*(dPhidx.*Deriv(:,1,Inod)+dPhidy.*Deriv(:,2,Inod)).*detJw;

        aw0term=dt*(1-theta)*aw0int.*fun(Inod).*detJw;
        aw1term=dt*theta*aw1int.*fun(Inod).*detJw;


        b1(:,Inod)=b1(:,Inod)+N0term+Vx0dN0dx+Vy0dN0dy+NTerm+PhiTerm+aw0term+aw1term ;

    end
end

% assemble right-hand side

R=sparseUA(neq,1);
for Inod=1:MUA.nod
    R=R+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b1(:,Inod),neq,1);
end


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

K=sparseUA(Iind,Jind,Xval,neq,neq);






end