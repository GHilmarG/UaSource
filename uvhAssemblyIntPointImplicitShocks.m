function   [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
    uvhAssemblyIntPointImplicitShocks(Iint,ndim,MUA,...
    bnod,hnod,unod,vnod,C,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dudt0nod,dvdt0nod,dadhnod,Bnod,Snod,rhonod,...
    CtrlVar,rhow,g,m,etaInt,exx,eyy,exy,Eint,Ronly,ca,sa,dt,...
    Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh)

% I've added here the rho terms in the mass-conservation equation
%

theta=CtrlVar.theta;


fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
%[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);

% tg3       241430610.225895
% theta     241430610.787565
% SUPG      241431244.153045
if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
else
    [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
end

%Deriv=MeshProp.Deriv{Iint} ; detJ=MeshProp.detJ{Iint};
% Deriv : MUA.Nele x dof x nod
%  detJ : MUA.Nele

dudt1nod=(unod-u0nod)/dt; dvdt1nod=(vnod-v0nod)/dt;

% values at integration this point
Hnod=Snod-Bnod;

hint=hnod*fun;
uint=unod*fun;
vint=vnod*fun;

if CtrlVar.CisElementBased
    Cint=C;
else
    Cint=C*fun;
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

hfint=rhow*Hint./rhoint;
kH=CtrlVar.kH;

Heint = HeavisideApprox(kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in consistent manner
HEint = HeavisideApprox(kH,hfint-hint,CtrlVar.Hh0);
deltaint=DiracDelta(kH,hint-hfint,CtrlVar.Hh0);      % i.e. deltaint must be the exact derivative of Heint
Deltaint=DiracDelta(kH,hfint-hint,CtrlVar.Hh0);      %  although delta is an even function...

%dintTest = HeavisideApprox(CtrlVar.kH,Sint-bint).*(Sint-bint);
%dintTest and dint are very similar
dint=HEint.*rhoint.*hint/rhow+Heint.*Hposint ;
%Dddhint=HEint.*rhoint/rhow+deltaint.*(Hint-rhoint.*hint/rhow);

Dddhint=HEint.*rhoint/rhow-Deltaint.*hint.*rhoint/rhow+deltaint.*Hposint;


[beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,m,Heint,CtrlVar);
etaint=etaInt(:,Iint) ;  % I could consider calculating this here


dhdx=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
dHdx=zeros(MUA.Nele,1); dHdy=zeros(MUA.Nele,1);
dBdx=zeros(MUA.Nele,1); dBdy=zeros(MUA.Nele,1);
dh0dx=zeros(MUA.Nele,1); dh0dy=zeros(MUA.Nele,1);
du0dx=zeros(MUA.Nele,1); du0dy=zeros(MUA.Nele,1);
dv0dx=zeros(MUA.Nele,1); dv0dy=zeros(MUA.Nele,1);
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
    
    du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
    du0dy=du0dy+Deriv(:,2,Inod).*u0nod(:,Inod);
    
    dv0dx=dv0dx+Deriv(:,1,Inod).*v0nod(:,Inod);
    dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
    
    dbdx=dbdx+Deriv(:,1,Inod).*bnod(:,Inod);
    dbdy=dbdy+Deriv(:,2,Inod).*bnod(:,Inod);
    
    drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
    drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);
    
end

qx1dx=rhoint.*exx(:,Iint).*hint+rhoint.*uint.*dhdx+drhodx.*uint.*hint;
qy1dy=rhoint.*eyy(:,Iint).*hint+rhoint.*vint.*dhdy+drhody.*vint.*hint;
qx0dx=rhoint.*du0dx.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*hint;
qy0dy=rhoint.*dv0dy.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*hint;

%%

l=2*sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));
beta=0.5;
tau0=beta*l./sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
%tau1=beta*l./sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);
%tau0=0; 
tau1=0;  % I have not included all corresponding derivaties in the left-hand side

detJw=detJ*MUA.weights(Iint);


for Inod=1:MUA.nod
    if ~Ronly
        for Jnod=1:MUA.nod
            
            
           
            
            
            Deu=Eint(:,Iint).*((2*exx(:,Iint)+eyy(:,Iint)).*Deriv(:,1,Jnod)+exy(:,Iint).*Deriv(:,2,Jnod));
            Dev=Eint(:,Iint).*((2*eyy(:,Iint)+exx(:,Iint)).*Deriv(:,2,Jnod)+exy(:,Iint).*Deriv(:,1,Jnod));
            % E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p
            E11=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Deu.*Deriv(:,1,Inod)+2*hint.*exy(:,Iint).*Deu.*Deriv(:,2,Inod);
            E12=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Dev.*Deriv(:,1,Inod)+2*hint.*exy(:,Iint).*Dev.*Deriv(:,2,Inod);
            E22=  hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Dev.*Deriv(:,2,Inod)+2*hint.*exy(:,Iint).*Dev.*Deriv(:,1,Inod);
            E21=  hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Deu.*Deriv(:,2,Inod)+2*hint.*exy(:,Iint).*Deu.*Deriv(:,1,Inod);
            
            
            Kxu(:,Inod,Jnod)=Kxu(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +beta2int.*fun(Jnod).*fun(Inod)...
                +E11...
                +Dbeta2Duuint.*fun(Jnod).*fun(Inod)).*detJw;
            
            Kyv(:,Inod,Jnod)=Kyv(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +beta2int.*fun(Jnod).*fun(Inod)...
                +E22...
                +Dbeta2Dvvint.*fun(Jnod).*fun(Inod)).*detJw ;
            
            Kxv(:,Inod,Jnod)=Kxv(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                +E12...
                +Dbeta2Duvint.*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
            
            
            Kyu(:,Inod,Jnod)=Kyu(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                +E21...
                +Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
            
            %   t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdx.*fun(Inod)+ rhoint.*g.*hint.*sa.*fun(Inod);
            %  t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
            %  Dddhint=HEint.*rhoint/rhow+deltaint.*Hposint-Deltaint.*hint.*rhoint/rhow;
            
            Kxh(:,Inod,Jnod)=Kxh(:,Inod,Jnod)...
                +(etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod).*fun(Jnod)...
                +etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod).*fun(Jnod)...
                +deltaint.*beta2int.*uint.*fun(Inod).*fun(Jnod)...
                +ca*g*rhoint.*Heint.*dBdx.*fun(Inod).*fun(Jnod)...                           % t1
                +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdx.*fun(Inod).*fun(Jnod)... ; % t1
                -sa*g*rhoint.*fun(Inod).*fun(Jnod)...                                        % t1
                -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,1,Inod).*fun(Jnod)...  ;    % t2
                ).*detJw;
            
            %-ca*g*rhoint.*(hint-dint.*rhoint.*HEint/rhow).*Deriv(:,1,Inod).*fun(Jnod)).*detJw;
            
            Kyh(:,Inod,Jnod)=Kyh(:,Inod,Jnod)...
                +(etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod).*fun(Jnod)...
                +etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod).*fun(Jnod)...
                +deltaint.*beta2int.*vint.*fun(Inod).*fun(Jnod)...
                +ca*g*rhoint.*Heint.*dBdy.*fun(Inod).*fun(Jnod)...                           % t1
                +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdy.*fun(Inod).*fun(Jnod)... ; % t1
                -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,2,Inod).*fun(Jnod)...  ;    % t2
                ).*detJw;
            
            %             SUPG=fun(Inod)+theta.*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
            %                 +(1-theta).*tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
            %
            N=(dh0dx.*dh0dx+dh0dy.*dh0dy);
            temp=(u0int.*dh0dx+v0int.*dh0dy)./N;
            temp(N<1e-5)=0;
            
            Shock=fun(Inod)+      theta.*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
                +(1-theta).* tau0.*temp.*(dh0dx.*Deriv(:,1,Inod)+dh0dy.*Deriv(:,2,Inod));
            
            
            
            Khu(:,Inod,Jnod)=Khu(:,Inod,Jnod)...
                +theta*(rhoint.*dhdx.*fun(Jnod)+drhodx.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,1,Jnod))...
                .*Shock.*detJw;
            
            
            Khv(:,Inod,Jnod)=Khv(:,Inod,Jnod)...
                +theta*(rhoint.*dhdy.*fun(Jnod)+drhody.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,2,Jnod))...
                .*Shock.*detJw;
            
            Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                +(rhoint.*fun(Jnod)/dt...
                -theta*rhoint.*dadhint.*fun(Jnod)...
                +theta*rhoint.*fun(Jnod).*h1barr/lambda_h...
                +theta.*(rhoint.*exx(:,Iint).*fun(Jnod)+drhodx.*uint.*fun(Jnod)+rhoint.*uint.*Deriv(:,1,Jnod)+...
                         rhoint.*eyy(:,Iint).*fun(Jnod)+drhody.*vint.*fun(Jnod)+rhoint.*vint.*Deriv(:,2,Jnod)))...
                .*Shock.*detJw;
            
        end
    end
    
    % note R=T-F;
    %  dR/dh  dh = -R
    %  dT/dh-dF/dh=-T+F  or dF/dh-dT/dh=T-F
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdx.*fun(Inod)+ rhoint.*g.*hint.*sa.*fun(Inod);
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
    t3=hint.*etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod);
    t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod);
    t5=beta2int.*uint.*fun(Inod);
    
    Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
    Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdy.*fun(Inod);
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,2,Inod);
    t3=hint.*etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod);
    t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod);
    t5=beta2int.*vint.*fun(Inod);
    
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
    
    %SUPG=fun(Inod)+      theta.*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
    %    +(1-theta).* tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
    
    % (u dhdx + v dhdy ) (dhdx ,dhdy)/ |grad h|^2
    % u0int*dh-dxi
    %
    
    
    N=(dh0dx.*dh0dx+dh0dy.*dh0dy);
    temp=(u0int.*dh0dx+v0int.*dh0dy)./N;
    
    temp(N<1e-5)=0;
    
    
    Shock=fun(Inod)+      theta.*tau1.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod))...
        +(1-theta).* tau0.*temp.*(dh0dx.*Deriv(:,1,Inod)+dh0dy.*Deriv(:,2,Inod));
        
    diffusion=(dh0dx.*Deriv(:,1,Inod)+dh0dy.*Deriv(:,2,Inod));  % does diffuse nicely
    
    qterm= (theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy).*Shock;
    
    D=1e7*l.*(1-Heint);
    %D=1e10*(1-Heint);
    
    qterm=qterm+D.*diffusion;
    
    dhdt=  rhoint.*((h0int-hint)/dt+(1-theta)*h0barr+theta*h1barr).*Shock;
    accterm=  rhoint.*((1-theta)*a0int+theta*a1int).*Shock;
    th=-dhdt;
    fh=  accterm - qterm;
    
    % R is calculated as R=th-fh  and then I solve K x = -R
    % thus: th has opposite sign but fh not

    Th(:,Inod)=Th(:,Inod)+th.*detJw;
    Fh(:,Inod)=Fh(:,Inod)+fh.*detJw;
    
    
end

end
