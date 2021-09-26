
function [UserVar,f0,K,dFdt]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,h0,rho,ub0,vb0,as0,ab0,h1,ub1,vb1,as1,ab1,das1dh,dab1dh)

% Assembly
%
%   K dh =-f0 
%
% dFdt is the matrix F in d deltah/dt = F deltah
% This matrix can be used to assess (linear) stability, from eigenvalues of M\dFdt
% 

narginchk(16,16)
nargoutchk(2,4)

nOut=nargout;

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.hTheta;
dt=CtrlVar.dt;




a1=as1+ab1;
a0=as0+ab0;
da1dh=das1dh+dab1dh;


h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(h1(MUA.connectivity,1),MUA.Nele,MUA.nod);

a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);

da1dhnod=reshape(da1dh(MUA.connectivity,1),MUA.Nele,MUA.nod);


ub0nod=reshape(ub0(MUA.connectivity,1),MUA.Nele,MUA.nod);
ub1nod=reshape(ub1(MUA.connectivity,1),MUA.Nele,MUA.nod);

vb0nod=reshape(vb0(MUA.connectivity,1),MUA.Nele,MUA.nod);
vb1nod=reshape(vb1(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


Khh=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFdt=zeros(MUA.Nele,MUA.nod,MUA.nod);
Rh=zeros(MUA.Nele,MUA.nod);


l=sqrt(2*MUA.EleAreas);
% vector over all elements for each  integration point

for Iint=1:MUA.nip  %Integration points
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    h0int=h0nod*fun;
    h1int=h1nod*fun;
    
    ub0int=ub0nod*fun; vb0int=vb0nod*fun;
    ub1int=ub1nod*fun; vb1int=vb1nod*fun;
    a0int=a0nod*fun;
    a1int=a1nod*fun;
    da1dhint=da1dhnod*fun;
    
    rhoint=rhonod*fun;
    
    % derivatives at one integration point for all elements
    Deriv1=squeeze(Deriv(:,1,:)) ;
    Deriv2=squeeze(Deriv(:,2,:)) ;
    
    
    exx0=zeros(MUA.Nele,1);
    eyy0=zeros(MUA.Nele,1);
    
    exx1=zeros(MUA.Nele,1);
    eyy1=zeros(MUA.Nele,1);
    
    drhodx=zeros(MUA.Nele,1); drhody=zeros(MUA.Nele,1);
    dh1dx=zeros(MUA.Nele,1); dh1dy=zeros(MUA.Nele,1);
    dh0dx=zeros(MUA.Nele,1); dh0dy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dh1dx=dh1dx+Deriv1(:,Inod).*h1nod(:,Inod);
        dh1dy=dh1dy+Deriv2(:,Inod).*h1nod(:,Inod);
        dh0dx=dh0dx+Deriv1(:,Inod).*h0nod(:,Inod);
        dh0dy=dh0dy+Deriv2(:,Inod).*h0nod(:,Inod);
        
        exx0=exx0+Deriv1(:,Inod).*ub0nod(:,Inod);
        eyy0=eyy0+Deriv2(:,Inod).*vb0nod(:,Inod);
        
        
        drhodx=drhodx+Deriv1(:,Inod).*rhonod(:,Inod);
        drhody=drhody+Deriv2(:,Inod).*rhonod(:,Inod);
        
        exx1=exx1+Deriv1(:,Inod).*ub1nod(:,Inod);
        eyy1=eyy1+Deriv2(:,Inod).*vb1nod(:,Inod);
        
        
        
    end
    
    
    detJw=detJ*MUA.weights(Iint);
    
    
    speed0=sqrt(ub0int.*ub0int+vb0int.*vb0int+CtrlVar.SpeedZero^2);
    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.h.SUPG.tau) ;
    tauSUPGint=CtrlVar.SUPG.beta0*tau;



    q1xdx=rhoint.*exx1.*h1int+rhoint.*ub1int.*dh1dx+drhodx.*ub1int.*h1int;
    q1ydy=rhoint.*eyy1.*h1int+rhoint.*vb1int.*dh1dy+drhody.*vb1int.*h1int;
    q0xdx=rhoint.*exx0.*h0int+rhoint.*ub0int.*dh0dx+drhodx.*ub0int.*h0int;
    q0ydy=rhoint.*eyy0.*h0int+rhoint.*vb0int.*dh0dy+drhody.*vb0int.*h0int;



    for Inod=1:MUA.nod


        SUPG=fun(Inod)+CtrlVar.h.SUPG.Use*tauSUPGint.*(ub0int.*Deriv1(:,Inod)+vb0int.*Deriv2(:,Inod));



        if nOut>2
            for Jnod=1:MUA.nod

                Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                    +(rhoint.*fun(Jnod)...
                    -dt*theta*rhoint.*da1dhint.*fun(Jnod)...
                    +dt*theta.*(rhoint.*exx1.*fun(Jnod)+drhodx.*ub1int.*fun(Jnod)+rhoint.*ub1int.*Deriv1(:,Jnod)...
                    +rhoint.*eyy1.*fun(Jnod)+drhody.*vb1int.*fun(Jnod)+rhoint.*vb1int.*Deriv2(:,Jnod)))...
                    .*SUPG.*detJw;

            end
        end


        if nOut>3
            for Jnod=1:MUA.nod

                dFdt(:,Inod,Jnod)=dFdt(:,Inod,Jnod)...
                    +(...
                    +theta*rhoint.*da1dhint.*fun(Jnod)...
                    -theta.*(rhoint.*exx1.*fun(Jnod)+drhodx.*ub1int.*fun(Jnod)+rhoint.*ub1int.*Deriv1(:,Jnod)...
                    -rhoint.*eyy1.*fun(Jnod)+drhody.*vb1int.*fun(Jnod)+rhoint.*vb1int.*Deriv2(:,Jnod)))...
                    .*SUPG.*detJw./rhoint;


            end
        end




        % Note, I solve: LSH  \phi  = - RHS
        qterm=  dt*(theta*q1xdx+(1-theta)*q0xdx+theta*q1ydy+(1-theta)*q0ydy);
        dhdt=  rhoint.*(h1int-h0int);
        accterm=  -dt*rhoint.*((1-theta)*a0int+theta*a1int);

        rh= dhdt + qterm + accterm ;

        % I solve K h = -Rh

        %
        Rh(:,Inod)=Rh(:,Inod)+rh.*SUPG.*detJw;


    end
end
%% assemble right-hand side

f0=sparseUA(neq,1);

for Inod=1:MUA.nod
    f0=f0+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Rh(:,Inod),neq,1);
end
%%




if nargout>2

    Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
    Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
    Kval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

    istak=0;
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod
            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
            Kval(istak+1:istak+MUA.Nele)=Khh(:,Inod,Jnod);
            istak=istak+MUA.Nele;
        end
    end

    K=sparseUA(Iind,Jind,Kval,neq,neq);

    if nargin>3
        dFdtVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                dFdtVal(istak+1:istak+MUA.Nele)=dFdt(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
        end
        dFdt=sparseUA(Iind,Jind,dFdtVal,neq,neq);
    end

end



end