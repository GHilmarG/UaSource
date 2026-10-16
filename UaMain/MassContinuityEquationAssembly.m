






function [UserVar,f0,K,dFdt]=MassContinuityEquationAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,BCs1)




%%
%
% [UserVar,f0,K,dFdt]=MassContinuityEquationAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)
%
% $$\rho \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \, \mathbf{v} h )  =  \rho \, a(h)$$
%
% $$a(h)$$ is a function of $h$ when using the level set with automated mass-balance feedback, and if using the thickness
%  barrier term to (approximately) enforce positive ice thicknesses.
%
% Assembly
%
%   K dh =-f0
%
% dFdt is the matrix F in d deltah/dt = F deltah This matrix can be used to assess (linear) stability, from eigenvalues of
% M\dFdt
%
%%


narginchk(6,7)
nargoutchk(2,4)

nOut=nargout;

% Nodes with thickness boundary conditions or active thickness constraints (8 Oct 2026). These nodes are excluded from the
% thickness penalty: with the integration-point penalty, elements containing such nodes are excluded, and with the nodal penalty
% (CtrlVar.ThicknessPenaltyNodal=true) the nodes themselves. If BCs1 is not given, no nodes are excluded.
if nargin<7 || isempty(BCs1)
    hConstrainedNodes=[];
else
    hConstrainedNodes=[BCs1.hFixedNode(:);BCs1.hPosNode(:)];
end
PenaltyNodal=isfield(CtrlVar,"ThicknessPenaltyNodal") && CtrlVar.ThicknessPenaltyNodal ;
PenaltyElementMask=~any(ismember(MUA.connectivity,hConstrainedNodes),2) ;   % Nele x 1, false for elements containing an h-constrained node

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.theta;
dt=CtrlVar.dt;


%% Mass-balance geometry feedback (8 Oct 2026), done in the same way as in uvhMatrixAssembly.
% For CtrlVar.MassBalanceGeometryFeedback=3 the mass balance is updated here using the current thickness (F1), at the end of the
% time step, and the derivative of the mass balance with respect to thickness is included in K (and dFdt). The derivative is
% treated as a nodal quantity: a_int=sum_k N_k a_k(h_k), so that d a_int/d h_j = N_j da_j/dh_j (see dadhnod below).
% For CtrlVar.MassBalanceGeometryFeedback=0 the mass balance is used as given, and no derivative is included.
% If the mass balance is evaluated at integration points ("-int-"), it is evaluated there directly as a function of the
% integration-point thickness, and the derivative is then an integration-point quantity (da1dhint) and dadhnod is zero.
MassBalanceAtIntPoints=contains(CtrlVar.MassBalance.Evaluation,"-int-");
if CtrlVar.MassBalanceGeometryFeedback==3 && ~MassBalanceAtIntPoints
    CtrlVar.time=CtrlVar.time+CtrlVar.dt;
    [UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);
    CtrlVar.time=CtrlVar.time-CtrlVar.dt;
    dadh=F1.dasdh+F1.dabdh;
else
    dadh=zeros(MUA.Nnodes,1);
end

a1=F1.as+F1.ab;
a0=F0.as+F0.ab;


h0nod=reshape(F0.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(F1.h(MUA.connectivity,1),MUA.Nele,MUA.nod);

a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);

dadhnod=reshape(dadh(MUA.connectivity,1),MUA.Nele,MUA.nod);   % nodal derivative of the (user-defined) mass balance


ub0nod=reshape(F0.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
ub1nod=reshape(F1.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);

vb0nod=reshape(F0.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
vb1nod=reshape(F1.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(F1.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

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

if CtrlVar.LevelSetMethod  &&  CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback

    if isempty(F1.LSFMask)
        F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
    end

    LSFMask=F1.LSFMask.NodesOut ; % This is the 'strictly' definition
    LSFMasknod=reshape(LSFMask(MUA.connectivity,1),MUA.Nele,MUA.nod);
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
        F1int.ab=nan(size(F1int.h)) ;

        GF0nodInt=GF0node*fun;
        GF1nodInt=GF1node*fun;

        F0int.GF.node=GF0nodInt;
        F1int.GF.node=GF1nodInt;

        MUAint.Nnodes=numel(h0int);

        if nargout("DefineMassBalance")==5
            [UserVar,as1int,ab1int,dasdhint,dabdhint]=DefineMassBalance(UserVar,CtrlVar,MUAint,F1int) ;
        else
            [UserVar,as1int,ab1int]=DefineMassBalance(UserVar,CtrlVar,MUAint,F1int) ;
            dasdhint=zeros(size(h1int)) ; dabdhint=zeros(size(h1int)) ;
        end
        if CtrlVar.MassBalanceGeometryFeedback~=3
            dasdhint=zeros(size(h1int)) ; dabdhint=zeros(size(h1int)) ;
        end
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
        da1dhint=zeros(MUA.Nele,1);   % the derivative of the user-defined mass balance is treated as a nodal quantity (dadhnod)

    end

    if CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback


        LM=LSFMasknod*fun;
        [abLSF,dadhLSF]=LevelSetMethodMassBalanceFeedback(CtrlVar,LM,h1int) ;
        a1int=a1int+abLSF;
        da1dhint=da1dhint+dadhLSF ;

    else

        LM=false ; % Level set mask for melt not applied

    end

    if isfield(CtrlVar,"ThicknessPenalty")  && CtrlVar.ThicknessPenalty && ~PenaltyNodal   % integration-point penalty (nodal penalty: see end of routine)


        %%  New simpler implementation of a thickness barrier.
        % Similar to the implementation of the LevelSetMethodAutomaticallyApplyMassBalanceFeedback the idea here is to directly
        % modify the mass-balance, a, and the da/dh rather than adding in new separate terms to the mass balance equation

        [aPenalty1,daPenaltydh1]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,h1int) ;
        aPenalty1=aPenalty1.*PenaltyElementMask ;          % no penalty in elements containing h-constrained nodes (8 Oct 2026)
        daPenaltydh1=daPenaltydh1.*PenaltyElementMask ;
        a1int=a1int+aPenalty1;
        da1dhint=da1dhint+daPenaltydh1 ;


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
    tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.uvh.SUPG.tau,CtrlVar.uvh.SUPG.tauMultiplier) ;
    tauSUPGint=CtrlVar.SUPG.beta0*tau;



    q1xdx=rhoint.*exx1.*h1int+rhoint.*ub1int.*dh1dx+drhodx.*ub1int.*h1int;
    q1ydy=rhoint.*eyy1.*h1int+rhoint.*vb1int.*dh1dy+drhody.*vb1int.*h1int;
    q0xdx=rhoint.*exx0.*h0int+rhoint.*ub0int.*dh0dx+drhodx.*ub0int.*h0int;
    q0ydy=rhoint.*eyy0.*h0int+rhoint.*vb0int.*dh0dy+drhody.*vb0int.*h0int;



    for Inod=1:MUA.nod


        SUPG=fun(Inod)+0.5*tauSUPGint.*(ub0int.*Deriv1(:,Inod)+vb0int.*Deriv2(:,Inod));   % same SUPG weighting as in the uvh assembly (uvhNodalLoopSSTREAM_v2)



        if nOut>2
            for Jnod=1:MUA.nod

                Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                    +(rhoint.*fun(Jnod)...
                    -dt*theta*rhoint.*(da1dhint+dadhnod(:,Jnod)).*fun(Jnod)...
                    +dt*theta.*(rhoint.*exx1.*fun(Jnod)+drhodx.*ub1int.*fun(Jnod)+rhoint.*ub1int.*Deriv1(:,Jnod)...
                               +rhoint.*eyy1.*fun(Jnod)+drhody.*vb1int.*fun(Jnod)+rhoint.*vb1int.*Deriv2(:,Jnod)))...
                    .*SUPG.*detJw;

            end
        end


        if nOut>3
            for Jnod=1:MUA.nod

                dFdt(:,Inod,Jnod)=dFdt(:,Inod,Jnod)...
                    +(...
                    +rhoint.*(da1dhint+dadhnod(:,Jnod)).*fun(Jnod)...  % no theta for continuous time derivative, theta is an averaging ratio of two discrete time steps, 
                    -...                            % but does not appear in the continuous formulation of the equations 
                    ( ...
                     rhoint.*exx1.*fun(Jnod)+drhodx.*ub1int.*fun(Jnod)+rhoint.*ub1int.*Deriv1(:,Jnod)...
                    +rhoint.*eyy1.*fun(Jnod)+drhody.*vb1int.*fun(Jnod)+rhoint.*vb1int.*Deriv2(:,Jnod) ...   % oops, there was a sign-error here, corrected on 8 Oct, 2026.
                    ) ...
                    )...
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
        accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int);

        rh= dhdt + qterm - accterm ;

        % I solve K h = -rh

        %
        Rh(:,Inod)=Rh(:,Inod)+rh.*SUPG.*detJw;


    end
end
%% assemble right-hand side

f0=sparse(neq,1);

for Inod=1:MUA.nod
    f0=f0+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Rh(:,Inod),neq,1);
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

    K=sparse(Iind,Jind,Kval,neq,neq);

    if nargout>3
        dFdtVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                dFdtVal(istak+1:istak+MUA.Nele)=dFdt(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
        end
        dFdt=sparse(Iind,Jind,dFdtVal,neq,neq);
    end

end


%% Lumped nodal thickness penalty (8 Oct 2026), done in the same way as in uvhMatrixAssembly (see comments there).
% The penalty is a nodal term that depends only on the nodal thickness, with HRZ-lumped nodal weights. Nodes with thickness
% boundary conditions or active thickness constraints are excluded. Added to f0, K and dFdt. As elsewhere in this routine, and as
% in the uvh assembly, CtrlVar.theta and the density F1.rho are used. In dFdt, which is the derivative of the continuous-time rate of
% change, there is no theta.
if isfield(CtrlVar,"ThicknessPenalty") && CtrlVar.ThicknessPenalty && PenaltyNodal

    MeDiag=zeros(MUA.Nele,MUA.nod) ; EleArea=zeros(MUA.Nele,1) ;
    for Iint=1:MUA.nip
        funI=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;    % nod x 1
        detJwI=MUA.DetJ(:,Iint)*MUA.weights(Iint) ;        % Nele x 1
        MeDiag=MeDiag+detJwI.*(funI.^2).' ;
        EleArea=EleArea+detJwI ;
    end
    mEle=MeDiag.*(EleArea./sum(MeDiag,2)) ;                % HRZ: element weights sum to element area
    mLump=accumarray(MUA.connectivity(:),mEle(:),[MUA.Nnodes 1]) ;

    wPen=ones(MUA.Nnodes,1) ; wPen(hConstrainedNodes)=0 ;
    [aPen,daPendh]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,F1.h) ;
    mw=mLump.*wPen ;
    f0=f0-sparse(dt*theta*F1.rho(:).*mw.*aPen) ;
    if nargout>2
        K=K+spdiags(-dt*theta*F1.rho(:).*mw.*daPendh,0,neq,neq) ;
    end
    if nargout>3
        dFdt=dFdt+spdiags(mw.*daPendh,0,neq,neq) ;
    end

end



end