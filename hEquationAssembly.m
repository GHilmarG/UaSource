
function [UserVar,kv,rh]=hEquationAssembly(UserVar,CtrlVar,MUA,u,v,a)

% Note: Assembly for the linear  h-equation:
%  d (u h)/dx + d (v h)/dy - div (kappa grad h) =  a 
%
% where h is the unknown.
%
%
%%

%%
%
% $$  \partial  ( u h)/ \partial x + \partial  ( v h)/\partial y - \nabla \cdot (\kappa \nabla h ) = a$$
%
%
% Isotropic:
%
%  $$ \kappa = k_{\mathrm{iso}} \; \mathbf{1}  $$ 
%
% Cross wind:
%
%  $$ \kappa = k_{\mathrm{cross}} \; ( \mathbf{1} - \hat{\mathbf{n}} \otimes \hat{\mathbf{n}})  $$ 
%
% Along wind:
%
%  $$ \kappa = k_{\mathrm{along}} \; (\hat{\mathbf{n}} \otimes \hat{\mathbf{n}} ) $$ 
%
% where $\hat{\mathbf{n}}$ is a unit vector pointing along the local flow direction.
% 
%%


narginchk(6,6)

ndim=2; dof=1; neq=dof*MUA.Nnodes;


kIso=CtrlVar.hEq.kIso;
kAlong=CtrlVar.hEq.kAlong;
kCross=CtrlVar.hEq.kCross;



unod=reshape(u(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
vnod=reshape(v(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
anod=reshape(a(MUA.connectivity,1),MUA.Nele,MUA.nod);


if isscalar(kIso)
    kIso=kIso+zeros(MUA.Nnodes,1);
end

if isscalar(kAlong)
    kAlong=kAlong+zeros(MUA.Nnodes,1);
end


if isscalar(kCross)
    kCross=kCross+zeros(MUA.Nnodes,1);
end





kIso=reshape(kIso(MUA.connectivity,1),MUA.Nele,MUA.nod);
kAlong=reshape(kAlong(MUA.connectivity,1),MUA.Nele,MUA.nod);
kCross=reshape(kCross(MUA.connectivity,1),MUA.Nele,MUA.nod);

d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
b1=zeros(MUA.Nele,MUA.nod);

l=sqrt(2*MUA.EleAreas);

% vector over all elements for each integration point
for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);

    % Deriv : Nele x dof x nod
    %  detJ : Nele

    % values at integration point


    uint=unod*fun;
    vint=vnod*fun;
    aint=anod*fun;

    kIsoint=kIso*fun;
    kAlongint=kAlong*fun;
    kCrossint=kCross*fun;


    speed=sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero) ; 
    nx=uint./speed; 
    ny=vint./speed;
    
    dudx=zeros(MUA.Nele,1); 
    dvdy=zeros(MUA.Nele,1); 
    
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);

    end
    
    detJw=detJ*MUA.weights(Iint);
        
    speed=sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);
    tau=SUPGtau(CtrlVar,speed,l,CtrlVar.dt,CtrlVar.Tracer.SUPG.tau) ;
    tauSUPGint=CtrlVar.SUPG.beta0*tau;

    for Inod=1:MUA.nod

        SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod));
        SUPGdetJw=SUPG.*detJw;

        for Jnod=1:MUA.nod

            hduxdvy=fun(Jnod).*(dudx+dvdy).*SUPGdetJw ;
            udhdxvdhdy=(uint.*Deriv(:,1,Jnod)+vint.*Deriv(:,2,Jnod)).*SUPGdetJw;


            % Isotropic diffusion term
            Diso=(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)).*detJw;

            % along-flow diffusion term
            Dalong=(nx.*Deriv(:,1,Jnod)+ny.*Deriv(:,2,Jnod)).*(nx.*Deriv(:,1,Inod)+ny.*Deriv(:,2,Inod)).*detJw;

            

            Diffusion=kIsoint.*Diso+kAlongint.*Dalong+kCrossint.*(Diso-Dalong) ;


            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+hduxdvy+udhdxvdhdy+Diffusion;

        end

        aterm=aint.*SUPGdetJw;
        b1(:,Inod)=b1(:,Inod)+aterm;

    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b1(:,Inod),neq,1);
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

kv=sparseUA(Iind,Jind,Xval,neq,neq);


end