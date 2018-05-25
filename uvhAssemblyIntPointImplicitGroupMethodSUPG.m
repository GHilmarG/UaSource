function   [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
    uvhAssemblyIntPointImplicitGroupMethodSUPG(Iint,ndim,MUA,...
    bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
    uonod,vonod,Conod,monod,uanod,vanod,Canod,manod,...
    CtrlVar,rhow,g,Ronly,ca,sa,dt,...
    Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh)

narginchk(52,52)


error('fsad')

% I've added here the rho terms in the mass-conservation equation
%
%
%  despite their names,  unod and Cnod can be either nodal or element variables
%
%

theta=CtrlVar.theta;

% SUPG weight
% fun(Inod)
% tau=l/ (2 |v|)


if ~CtrlVar.IncludeMelangeModelPhysics
    uoint=[];
    voint=[];
    Coint=[];
    moint=[];
    uaint=[];
    vaint=[];
    Caint=[];
    maint=[];
end

fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
else
    [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
end

%Deriv=MeshProp.Deriv{Iint} ; detJ=MeshProp.detJ{Iint};
% Deriv : MUA.Nele x dof x nod
%  detJ : MUA.Nele

% values at integration this point
Hnod=Snod-Bnod;

hint=hnod*fun;
uint=unod*fun;
vint=vnod*fun;


if CtrlVar.IncludeMelangeModelPhysics
    
    uoint=uonod*fun;
    voint=vonod*fun;
    
    uaint=uanod*fun;
    vaint=vanod*fun;
    
end

if CtrlVar.CisElementBased
    Cint=Cnod;
    mint=mnod;
    
    if CtrlVar.IncludeMelangeModelPhysics
        Coint=Conod;
        moint=monod;
        
        Caint=Canod;
        maint=manod;
    end
    
    
    
else
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    
    if CtrlVar.IncludeMelangeModelPhysics
        Coint=Conod*fun;
        moint=monod*fun;
        
        Caint=Canod*fun;
        maint=manod*fun;
    end
end


if CtrlVar.AGlenisElementBased
    AGlenint=AGlennod;
    nint=nnod;
else
    AGlenint=AGlennod*fun;
    AGlenint(AGlenint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    nint=nnod*fun;
end


h0int=h0nod*fun;
u0int=u0nod*fun;
v0int=v0nod*fun;

as0int=as0nod*fun;
ab0int=ab0nod*fun;
as1int=as1nod*fun;
ab1int=ab1nod*fun;
a1int=as1int+ab1int;
a0int=as0int+ab0int;
dadhint=dadhnod*fun;

Bint=Bnod*fun;
Sint=Snod*fun;
rhoint=rhonod*fun;



Hint=Sint-Bint;

if CtrlVar.ThicknessBarrier
    
    % using ThicknessBarrier I add fictitious accumulation term:
    %
    % gamma exp(-(h-h0)/l)
    % where:
    %       gamma=CtrlVar.ThicknessBarrierAccumulation
    %       l=CtrlVar.ThicknessBarrierThicknessScale
    %       h0=CtrlVar.ThickMin*CtrlVar.ThicknessBarrierMinThickMultiplier
    %
    lambda_h=CtrlVar.ThicknessBarrierThicknessScale;
    gamma_h=CtrlVar.ThicknessBarrierAccumulation;
    
    ThickBarrierMin=CtrlVar.ThickMin*CtrlVar.ThicknessBarrierMinThickMultiplier;
    
    argmax=log(realmax)/2;
    
    %       don't add (fictitious/barrier) mass term to last time step, only to the new one
    %        hTest=h0int ;
    %        arg0=-(hTest-MeshSizeMin)/lambda_h;
    %        arg0(arg0>argmax)=argmax;
    %        h0barr=gamma_h*exp(arg0)/lambda_h;
    h0barr=0;
    
    
    arg1=-(hint-ThickBarrierMin)/lambda_h;
    arg1(arg1>argmax)=argmax;
    h1barr=gamma_h*exp(arg1)/lambda_h;
    %
    
    %fprintf(' max(h1barr)=%-g max(h0barr)=%-g \n ',max(max(h1barr)),max(max(h0barr)))
else
    h1barr=0 ; h0barr=0; lambda_h=1;
end


Hposint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*Hint;


% calculate He and DiracDelta from integration point values of thickness

%% Group approach to floating mask
%
% 1) evaluate hfloat at nodes as : hf= rhow*(S-B)./rho
% 2) evaluate floating mask at nodes as : He=HeavisideApprox(h-hf)
% 3) evaluate spatial derivative of floating mask using usuall FE approach

hfint=rhow*Hint./rhoint;
kH=CtrlVar.kH;

Heint = HeavisideApprox(kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
HEint = HeavisideApprox(kH,hfint-hint,CtrlVar.Hh0);
deltaint=DiracDelta(kH,hint-hfint,CtrlVar.Hh0);      % i.e. deltaint must be the exact derivative of Heint
Deltaint=DiracDelta(kH,hfint-hint,CtrlVar.Hh0);      %  although delta is an even function...

dint=HEint.*rhoint.*hint/rhow+Heint.*Hposint ;
Dddhint=HEint.*rhoint/rhow-Deltaint.*hint.*rhoint/rhow+deltaint.*Hposint;



dhdx=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
dHdx=zeros(MUA.Nele,1); dHdy=zeros(MUA.Nele,1);
dBdx=zeros(MUA.Nele,1); dBdy=zeros(MUA.Nele,1);
dh0dx=zeros(MUA.Nele,1); dh0dy=zeros(MUA.Nele,1);



exx0=zeros(MUA.Nele,1);
eyy0=zeros(MUA.Nele,1);

exx=zeros(MUA.Nele,1);
eyy=zeros(MUA.Nele,1);
exy=zeros(MUA.Nele,1);

dbdx=zeros(MUA.Nele,1); dbdy=zeros(MUA.Nele,1);
drhodx=zeros(MUA.Nele,1); drhody=zeros(MUA.Nele,1);



% derivatives at integration points
for Inod=1:MUA.nod
    
    dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
    dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
    
    dHdx=dHdx+Deriv(:,1,Inod).*Hnod(:,Inod);
    dHdy=dHdy+Deriv(:,2,Inod).*Hnod(:,Inod);
    
    dBdx=dBdx+Deriv(:,1,Inod).*Bnod(:,Inod);
    dBdy=dBdy+Deriv(:,2,Inod).*Bnod(:,Inod);
    
    dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
    dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
    
    exx0=exx0+Deriv(:,1,Inod).*u0nod(:,Inod);  % exx0
    %du0dy=du0dy+Deriv(:,2,Inod).*u0nod(:,Inod);
    
    %dv0dx=dv0dx+Deriv(:,1,Inod).*v0nod(:,Inod);
    eyy0=eyy0+Deriv(:,2,Inod).*v0nod(:,Inod);
    
    dbdx=dbdx+Deriv(:,1,Inod).*bnod(:,Inod);
    dbdy=dbdy+Deriv(:,2,Inod).*bnod(:,Inod);
    
    drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
    drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);
    
    exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
    eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
    exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod));
    
    
end


[etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);

%uoint=[];voint=[];Coint=[] ;moint=[] ;uaint=[] ;vaint=[] ;Caint=[]; maint=[];
[taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh] = BasalDrag(CtrlVar,Heint,deltaint,hint,Bint,Hint,rhoint,rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint);



CtrlVar.GroupRepresentation=0;


qx1dx=rhoint.*exx.*hint+rhoint.*uint.*dhdx+drhodx.*uint.*hint;
qy1dy=rhoint.*eyy.*hint+rhoint.*vint.*dhdy+drhody.*vint.*hint;
qx0dx=rhoint.*exx0.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*hint;
qy0dy=rhoint.*eyy0.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*hint;





%% u=1 ; dt =1 ; l=1 ; tau=1/(u/l+l/(u*dt^2))

l=2*sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));

speed0=sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
speed1=sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);

%
% The SUPG adds a term to the weighting function on the form
%
%  N -> N+ N'
%
%  where
%
%   N'=tau \bm{u} \cdot \grad N
%     = tau (u dNdx + v dNdy)
%
% where tau is a parameter having the dimension time. In a transient run a
% possible choice for tau is simply dt However we would like the aditional SUPG
% term to go to zero as u->0 and as h->0, in fact we expect the `element Courant
% number' ECN defined as
%
%               ECN = dt |\bm{u}|/L,
%
% where u is some typical velocity and L a scale for element size to be of relevance.
%
%  Typical suggestions in the litterature are on the form
%
%  N'= kappa (L/2) u/|u| dNdx
%
% for a 1D situation on a regular grid, with]%
%
% kappa=cosh(Pe) -1/Pe
%
% and Pe=U L / (2 k)  , where k is the diffusivity constant.
%
% Huges et al suggest
%
%          tau= \frac{L}{2 |\bm{u}|}    (coth ( Pe) - 1/Pe )
%
% with Pe=|bm{u}| L / 2 k,  where k is the diffusion coefficient.
% It is unclear what the diffusion coefficient will be. In the hyperbolic
% limit of a k->0, kappa=cosh(Pe)-1/Pe -> 1 and the SUPG becomes
%
%  N'=   (L/2|u|)  \bm{u} \cdot \grad N
%
% This terms goes to zero with decreasing element size.
%
%
% Heuristic arguments
% suggest \alpha=ECN as this gives plausible limites, i.e. 1 for ECN->\infty and
% 0 for ECN-> 0.
%
%
% (L/speed)* coth( speed dt/L)  * (u dN/dx + v dN/dy)
% _
% Note: ((coth(x)-1/x)/x -> 1/3 for x ->0
%       ((coth(1/x)-x) x -> 1/3 for x -> infty
%       ((coth(x)-1/x)  -> 0 for x ->0
%       ((coth(x)-1/x)  -> 1 for x -> infty
%
% this term is zero for L->0, or speed->0, or dt-> 0
%
%       dN/dX for u-> infty ,  dt-> infty
%  1/3  dN/dX for L -> infty
%
ECN=100*speed0.*dt./l+eps; % This is the `Element Courant Number' ECN
kappa=coth(ECN)-1./ECN;
tau0=CtrlVar.SUPG.beta0*kappa.*l./speed0 ; % sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
tau1=CtrlVar.SUPG.beta1*kappa.*l./speed1 ; % sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);

f=rhoint.*(hint-h0int-dt*(1-theta)*h0barr-dt*theta*h1barr)+dt*theta*qx1dx+dt*(1-theta)*qx0dx+dt*theta*qy1dy+dt*(1-theta)*qy0dy-dt*rhoint.*((1-theta)*a0int+theta*a1int);
SUPGu=theta*CtrlVar.SUPG.beta1*l.*f.*(1./sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)-uint.*uint./(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2).^(3/2));
SUPGv=theta*CtrlVar.SUPG.beta1*l.*f.*(1./sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)-vint.*vint./(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2).^(3/2));


detJw=detJ*MUA.weights(Iint);


for Inod=1:MUA.nod
    if ~Ronly
        for Jnod=1:MUA.nod
            
            Deu=Eint.*((2*exx+eyy).*Deriv(:,1,Jnod)+exy.*Deriv(:,2,Jnod));
            Dev=Eint.*((2*eyy+exx).*Deriv(:,2,Jnod)+exy.*Deriv(:,1,Jnod));
            % E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p
            E11=  hint.*(4.*exx+2.*eyy).*Deu.*Deriv(:,1,Inod)+2*hint.*exy.*Deu.*Deriv(:,2,Inod);
            E12=  hint.*(4.*exx+2.*eyy).*Dev.*Deriv(:,1,Inod)+2*hint.*exy.*Dev.*Deriv(:,2,Inod);
            E22=  hint.*(4.*eyy+2.*exx).*Dev.*Deriv(:,2,Inod)+2*hint.*exy.*Dev.*Deriv(:,1,Inod);
            E21=  hint.*(4.*eyy+2.*exx).*Deu.*Deriv(:,2,Inod)+2*hint.*exy.*Deu.*Deriv(:,1,Inod);
            
            
            Kxu(:,Inod,Jnod)=Kxu(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +E11...
                +dtauxdu.*fun(Jnod).*fun(Inod)...   % +beta2int.*fun(Jnod).*fun(Inod)+Dbeta2Duuint.*fun(Jnod).*fun(Inod)...
                ).*detJw;
            
            Kyv(:,Inod,Jnod)=Kyv(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +E22...
                +dtauydv.*fun(Jnod).*fun(Inod)...   % +beta2int.*fun(Jnod).*fun(Inod)+Dbeta2Dvvint.*fun(Jnod).*fun(Inod)...
                ).*detJw ;
            
            Kxv(:,Inod,Jnod)=Kxv(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                +E12...
                +dtauxdv.*fun(Jnod).*fun(Inod)...   % +Dbeta2Duvint.*fun(Jnod).*fun(Inod)...    % beta derivative
                ).*detJw;
            
            
            Kyu(:,Inod,Jnod)=Kyu(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                +E21...
                +dtauydu.*fun(Jnod).*fun(Inod)...    % +Dbeta2Duvint*fun(Jnod).*fun(Inod)...
                ).*detJw;
            
            %   t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdx.*fun(Inod)+ rhoint.*g.*hint.*sa.*fun(Inod);
            %  t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
            %  Dddhint=HEint.*rhoint/rhow+deltaint.*Hposint-Deltaint.*hint.*rhoint/rhow;
            
            Kxh(:,Inod,Jnod)=Kxh(:,Inod,Jnod)...
                +(etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod).*fun(Jnod)...
                +etaint.*2.*exy.*Deriv(:,2,Inod).*fun(Jnod)...
                +dtauxdh.*fun(Inod).*fun(Jnod)... % +deltaint.*beta2int.*uint.*fun(Inod).*fun(Jnod)..
                +ca*g*rhoint.*Heint.*dBdx.*fun(Inod).*fun(Jnod)...                           % t1
                +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdx.*fun(Inod).*fun(Jnod)... ; % t1
                -sa*g*rhoint.*fun(Inod).*fun(Jnod)...                                        % t1
                -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,1,Inod).*fun(Jnod)...  ;    % t2
                ).*detJw;
            
            %-ca*g*rhoint.*(hint-dint.*rhoint.*HEint/rhow).*Deriv(:,1,Inod).*fun(Jnod)).*detJw;
            
            Kyh(:,Inod,Jnod)=Kyh(:,Inod,Jnod)...
                +(etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod).*fun(Jnod)...
                +etaint.*2.*exy.*Deriv(:,1,Inod).*fun(Jnod)...
                +dtauydh.*fun(Inod).*fun(Jnod)...   % +deltaint.*beta2int.*vint.*fun(Inod).*fun(Jnod)...
                +ca*g*rhoint.*Heint.*dBdy.*fun(Inod).*fun(Jnod)...                           % t1
                +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdy.*fun(Inod).*fun(Jnod)... ; % t1
                -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,2,Inod).*fun(Jnod)...  ;    % t2
                ).*detJw;
            
            SUPG=fun(Inod)+theta.*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
                +(1-theta).*tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
            
            
            dSUPGu=SUPGu.*fun(Jnod).*Deriv(:,1,Inod);
            dSUPGv=SUPGv.*fun(Jnod).*Deriv(:,2,Inod);
            
            Khu(:,Inod,Jnod)=Khu(:,Inod,Jnod)...
                +theta*(rhoint.*dhdx.*fun(Jnod)+drhodx.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,1,Jnod))...
                .*SUPG.*detJw*dt+dSUPGu.*detJw*dt;
            
            
            Khv(:,Inod,Jnod)=Khv(:,Inod,Jnod)...
                +theta*(rhoint.*dhdy.*fun(Jnod)+drhody.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,2,Jnod))...
                .*SUPG.*detJw*dt+dSUPGv.*detJw*dt;
            
            
            
            Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                +(rhoint.*fun(Jnod)...
                -dt*theta*rhoint.*dadhint.*fun(Jnod)...
                +dt*theta*rhoint.*fun(Jnod).*h1barr/lambda_h...
                +dt*theta.*(rhoint.*exx.*fun(Jnod)+drhodx.*uint.*fun(Jnod)+rhoint.*uint.*Deriv(:,1,Jnod)+...
                rhoint.*eyy.*fun(Jnod)+drhody.*vint.*fun(Jnod)+rhoint.*vint.*Deriv(:,2,Jnod)))...
                .*SUPG.*detJw;
            
        end
    end
    
    % note R=T-F;
    %  dR/dh  dh = -R
    %  dT/dh-dF/dh=-T+F  or dF/dh-dT/dh=T-F
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdx.*fun(Inod)+ rhoint.*g.*hint.*sa.*fun(Inod);
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
    t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);
    t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);
    t5=taux.*fun(Inod); % beta2int.*uint.*fun(Inod);
    
    Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
    Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdy.*fun(Inod);
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,2,Inod);
    t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
    t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
    t5=tauy.*fun(Inod); % beta2int.*vint.*fun(Inod);
    
    Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
    Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
    
    %qxx= d( uh)/dx
    
    % R=T-F
    
    % first-order Taylor terms
    % th=  (theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy).*fun(Inod);
    % fh=  ((h0int-hint)/dt+(1-theta)*(a0int+h0barr)+theta*(a1int+h1barr)).*fun(Inod);
    
    
    
    %         th1=  (theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy).*fun(Inod);
    %         th2=  ((h0int-hint)/dt+(1-theta)*h0barr+theta*h1barr).*fun(Inod);
    %         th=th1-th2;
    %         fh=  ((1-theta)*a0int+theta*a1int).*fun(Inod);
    
    % changed the def of th and fh (jan 2014)
    
    % fun(Inod) -> fun(Inod)+theta .* tau.*uint.*Deriv(:,2,Inod)+vint.*Deriv(:,2,Inod))
    %                    +(1-theta).* tau.*u0int.*Deriv(:,2,Inod)+v0int.*Deriv(:,2,Inod))
    %
    
    SUPG=fun(Inod)+    theta .*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
        +(1-theta).* tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
    
    
    qterm=  dt*(theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy).*SUPG;
    dhdt=  rhoint.*(h0int-hint+dt*(1-theta)*h0barr+dt*theta*h1barr).*SUPG;
    accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int).*SUPG;
    
    th=-dhdt;
    fh=  accterm - qterm;
    
    % R is calculated as R=th-fh  and then I solve K x = -R
    % thus: th has opposite sign but fh not
    % second and third-order Taylor terms
    
    %
    Th(:,Inod)=Th(:,Inod)+th.*detJw;
    Fh(:,Inod)=Fh(:,Inod)+fh.*detJw;
    
    
end

end
