
function   [UserVar,f0,K,dFdt]=MassContinuityEquationAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)


% [UserVar,f0,K,dFdt]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,h0,rho,ub0,vb0,as0,ab0,h1,ub1,vb1,as1,ab1,das1dh,dab1dh)

% Assembly
%
%   K dh =-f0
%
% dFdt is the matrix F in d deltah/dt = F deltah
% This matrix can be used to assess (linear) stability, from eigenvalues of M\dFdt
%

narginchk(6,6)
nargoutchk(2,4)

nOut=nargout;

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.hTheta;
dt=CtrlVar.dt;




a1=F1.as+F1.ab;
a0=F0.as+F0.ab;
da1dh=F1.dasdh+F1.dabdh;


h0nod=reshape(F0.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(F1.h(MUA.connectivity,1),MUA.Nele,MUA.nod);

a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);

da1dhnod=reshape(da1dh(MUA.connectivity,1),MUA.Nele,MUA.nod);


ub0nod=reshape(F0.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
ub1nod=reshape(F1.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);

vb0nod=reshape(F0.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
vb1nod=reshape(F1.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(F0.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

s0nod=reshape(F0.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
s1nod=reshape(F1.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

b0nod=reshape(F0.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
b1nod=reshape(F1.b(MUA.connectivity,1),MUA.Nele,MUA.nod);

coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);

GF0node=reshape(F0.GF.node(MUA.connectivity,1),MUA.Nele,MUA.nod);
GF1node=reshape(F1.GF.node(MUA.connectivity,1),MUA.Nele,MUA.nod);


Khh=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFdt=zeros(MUA.Nele,MUA.nod,MUA.nod);
Rh=zeros(MUA.Nele,MUA.nod);

if isempty(F1.dabdh)
    F1.dabdh=zeros(MUA.Nnodes,1) ;
end

if CtrlVar.LevelSetMethod  &&  CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback  && ~isempty(F1.LSF)

    if ~isempty(F1.LSF) 

        if isempty(F1.LSFMask)
            F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
        end

        LSFMask=F1.LSFMask.NodesOut ; % This is the 'strickly' definition
        LSFMasknod=reshape(LSFMask(MUA.connectivity,1),MUA.Nele,MUA.nod);


    end


end




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
    rhoint=rhonod*fun;

    

    if  contains(CtrlVar.MassBalance.Evaluation,"-int-")


        % This is not fully flexible as I have not implemented calculating the flotation mask at integration points
        xint=coox*fun;  % coordinates of this integration point for all elements
        yint=cooy*fun;

        F0int.x=xint;  F0int.y=yint;
        F1int.x=xint;  F1int.y=yint;

        s0int=s0nod*fun;
        s1int=s1nod*fun;

        b0int=b0nod*fun;
        b1int=b1nod*fun;


        F0int.h=h0int;
        F0int.s=s0int;
        F0int.b=b0int;
        F0int.rho=rhoint;
        F0int.rhow=F0.rhow;
        F0int.S=F0.S(1);

        F1int.h=h1int;
        F1int.s=s1int;
        F1int.b=b1int;
        F1int.rho=rhoint;
        F1int.rhow=F1.rhow;
        F1int.S=F1.S(1);

        F0int.as=nan(size(F0int.h)) ;
        F0int.ab=nan(size(F0int.h)) ;
        F1int.as=nan(size(F1int.h)) ;
        F1int.as=nan(size(F1int.h)) ;

        GF0nodInt=GF0node*fun;
        GF1nodInt=GF1node*fun;

        F0int.GF.node=GF0nodInt;
        F1int.GF.node=GF1nodInt;

        MUAint.Nnodes=numel(h0int);

        [UserVar,as1int,ab1int,dasdhint,dabdhint]=DefineMassBalance(UserVar,CtrlVar,MUAint,F1int) ;
        [UserVar,as0int,ab0int]=DefineMassBalance(UserVar,CtrlVar,MUAint,F0int) ;

        a0int=as0int+ab0int ;
        a1int=as1int+ab1int ;
        da1dhint=dasdhint+dabdhint;

        % a0intTest=a0nod*fun;
        % a1intTest=a1nod*fun;
        % da1dhintTest=da1dhnod*fun;
        %
        % [norm(a0int-a0intTest) norm(a1int-a1intTest) norm(da1dhintTest-da1dhint)]

    else

        a0int=a0nod*fun;
        a1int=a1nod*fun;
        da1dhint=da1dhnod*fun;

    end

    if CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback

        LM=LSFMasknod*fun;
        a1= CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffLin;
        a3= CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffCubic;

        hmin=CtrlVar.LevelSetMinIceThickness;

        abLSF =LM.* ( a1*(h1int-hmin)+a3*(h1int-hmin).^3) ;
        dadhLSF=LM.*(a1+3*a3*(h1int-hmin).^2) ;

        a1int=a1int+abLSF; 
        da1dhint=da1dhint+dadhLSF ;
    end
    



    % da1dhint=0;

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


                % dFdt = ( rho a - div ( rho v)   )/rho
                %
                %  A dh/dt = F
                %
                %   rho (h1 - h0)/dt  =   rho a - div (rho v)
                %
                %  < rho (h1 - h0)/dt  , \phi > =  < rho a - div (rho v)  , \phi >
                %
                %
                %

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