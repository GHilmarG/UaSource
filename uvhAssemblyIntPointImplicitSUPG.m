function   [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
    uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
    bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,qnod,muknod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
    Henod,deltanod,Hposnod,dnod,Dddhnod,...
    LSFMasknod,...
    uonod,vonod,Conod,monod,uanod,vanod,Canod,manod,...
    CtrlVar,rhow,g,Ronly,ca,sa,dt,...
    Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh)

narginchk(60,60)


% I've added here the rho terms in the mass-conservation equation
%
%  despite their names,  unod and Cnod can be either nodal or element variables
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


% if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
Deriv=MUA.Deriv(:,:,:,Iint);
detJ=MUA.DetJ(:,Iint);
% else
% [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
% end

%Deriv=MeshProp.Deriv{Iint} ; detJ=MeshProp.detJ{Iint};
% Deriv : MUA.Nele x dof x nod
%  detJ : MUA.Nele

% values at integration this point
% Hnod=Snod-Bnod;

hint=hnod*fun;
uint=unod*fun;
vint=vnod*fun;


if CtrlVar.IncludeMelangeModelPhysics
    
    uoint=uonod*fun;
    voint=vonod*fun;
    
    uaint=uanod*fun;
    vaint=vanod*fun;
    
end


Cint=Cnod*fun;
% Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
mint=mnod*fun;

if ~isempty(qnod)
    qint=qnod*fun;
else
    qint=[];
end


if ~isempty(muknod)
    mukint=muknod*fun;
else
    mukint=[];
end

if CtrlVar.IncludeMelangeModelPhysics
    Coint=Conod*fun;
    moint=monod*fun;
    
    Caint=Canod*fun;
    maint=manod*fun;
end

AGlenint=AGlennod*fun;
% AGlenint(AGlenint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
nint=nnod*fun;



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






if CtrlVar.LevelSetMethod  &&  CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback

    LM=LSFMasknod*fun;
    a1= CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffLin;
    a3= CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffCubic;
  
    hmin=CtrlVar.LevelSetMinIceThickness;

    abLSF =LM.* ( a1*(hint-hmin)+a3*(hint-hmin).^3) ;
    dadhLSF=LM.*(a1+3*a3*(hint-hmin).^2) ;

    a1int=a1int+abLSF; dadhint=dadhint+dadhLSF ;

else
    LM=false; % Level set mask for melt not applied
end

h1barr=0 ; h0barr=0; lambda_h=1;

if CtrlVar.ThicknessBarrier

    % using ThicknessBarrier I add fictitious accumulation term:
    %
    % gamma exp(-(h-h0)/l)
    % where:
    %       gamma=CtrlVar.ThicknessBarrierAccumulation
    %       l=CtrlVar.ThicknessBarrierThicknessScale
    %       h0=CtrlVar.ThickMin*CtrlVar.ThicknessBarrierMinThickMultiplier
    %
    %     lambda_h=CtrlVar.ThicknessBarrierThicknessScale;
    %     gamma_h=CtrlVar.ThicknessBarrierAccumulation;
    %
    %     ThickBarrierMin=CtrlVar.ThickMin*CtrlVar.ThicknessBarrierMinThickMultiplier;
    %
    %     argmax=log(realmax)/2;
    %     h0barr=0;
    %
    %
    %     arg1=-(hint-ThickBarrierMin)/lambda_h;
    %     arg1(arg1>argmax)=argmax;
    %     h1barr=gamma_h*exp(arg1)/lambda_h;


    %%  New simpler implementation of a thickness barrier.
    % Similar to the implementatoin of the LevelSetMethodAutomaticallyApplyMassBalanceFeedback
    % the idea here is to directly modify the mass-balance, a, and the da/dh rather than adding in new seperate terms to the mass
    % balance equation

    hmin=CtrlVar.ThickMin ;

    isThickTooSmall=hint<hmin ;

    % don't apply if already applied as a part of the level-set method
    isThickTooSmall=isThickTooSmall & ~LM ;
    a1= CtrlVar.ThicknessBarrierMassBalanceFeedbackCoeffLin;
    a3= CtrlVar.ThicknessBarrierMassBalanceFeedbackCoeffCubic;


    abThickMin =isThickTooSmall.* ( a1*(hint-hmin)+a3*(hint-hmin).^3) ;  % if thickness too small, then (hint-hmin) < 0, and ab > 0

    dadhThickMin=isThickTooSmall.*(a1+3*a3*(hint-hmin).^2) ;

    a1int=a1int+abThickMin; dadhint=dadhint+dadhThickMin ;

    %     nh=numel(find(isThickTooSmall)) ;
    %     if nh> 0
    %        fprintf("#%i \t ThickMin=%f \t max(abThickMin)=%f \n ",nh,min(hint(isThickTooSmall)),max(abThickMin))
    %     end
    %
end

Bint=Bnod*fun;
Sint=Snod*fun;
rhoint=rhonod*fun;
Hint=Sint-Bint;

hfint=rhow*Hint./rhoint;  % this is linear, so fine to evaluate at int in this manner


if CtrlVar.uvhGroupAssembly


    Heint=Henod*fun;  %
   % HEint=HEnod*fun;  %

    deltaint=deltanod*fun;
   % Deltaint=Deltanod*fun;

    Hposint=Hposnod*fun;

    dint=dnod*fun;        % dnod=HeavisideApprox(CtrlVar.kH,Hnod,CtrlVar.Hh0).*(Snod-bnod);  % draft
    Dddhint=Dddhnod*fun;  % Here I am interpolating the derivative calculated at nodes, to the int points


else

    % calculate He and DiracDelta from integration point values of thickness

   
 

    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
    HEint = HeavisideApprox(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);
    
    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);      % i.e. deltaint must be the exact derivative of Heint
    Deltaint=DiracDelta(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);      %  although delta is an even function...

    Hposint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*Hint;

    dint=HEint.*rhoint.*hint/rhow+Heint.*Hposint ;  % definition of d
    Dddhint=HEint.*rhoint/rhow-Deltaint.*hint.*rhoint/rhow+deltaint.*Hposint; % derivative of dint with respect to hint

end


dhdx=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
%dHdx=zeros(MUA.Nele,1); dHdy=zeros(MUA.Nele,1);
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

Deriv1=squeeze(Deriv(:,1,:)) ;
Deriv2=squeeze(Deriv(:,2,:)) ;

for Inod=1:MUA.nod
    
    dhdx=dhdx+Deriv1(:,Inod).*hnod(:,Inod);
    dhdy=dhdy+Deriv2(:,Inod).*hnod(:,Inod);
    
    %    dHdx=dHdx+Deriv1(:,Inod).*Hnod(:,Inod);
    %    dHdy=dHdy+Deriv2(:,Inod).*Hnod(:,Inod);
    
    dBdx=dBdx+Deriv1(:,Inod).*Bnod(:,Inod);
    dBdy=dBdy+Deriv2(:,Inod).*Bnod(:,Inod);
    
    dh0dx=dh0dx+Deriv1(:,Inod).*h0nod(:,Inod);
    dh0dy=dh0dy+Deriv2(:,Inod).*h0nod(:,Inod);
    
    exx0=exx0+Deriv1(:,Inod).*u0nod(:,Inod);  % exx0
    eyy0=eyy0+Deriv2(:,Inod).*v0nod(:,Inod);
    
    dbdx=dbdx+Deriv1(:,Inod).*bnod(:,Inod);
    dbdy=dbdy+Deriv2(:,Inod).*bnod(:,Inod);
    
    drhodx=drhodx+Deriv1(:,Inod).*rhonod(:,Inod);
    drhody=drhody+Deriv2(:,Inod).*rhonod(:,Inod);
    
    exx=exx+Deriv1(:,Inod).*unod(:,Inod);
    eyy=eyy+Deriv2(:,Inod).*vnod(:,Inod);
    exy=exy+0.5*(Deriv1(:,Inod).*vnod(:,Inod) + Deriv2(:,Inod).*unod(:,Inod));
    
    
end



[etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);

%uoint=[];voint=[];Coint=[] ;moint=[] ;uaint=[] ;vaint=[] ;Caint=[]; maint=[];
[taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh] = ...
    BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint,qint,g,mukint);



%% u=1 ; dt =1 ; l=1 ; tau=1/(u/l+l/(u*dt^2))

%l=2*sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));

l=sqrt(2*MUA.EleAreas);  % this I could take outside the int loop

speed0=sqrt(u0int.*u0int+v0int.*v0int+CtrlVar.SpeedZero^2);
% speed1=sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2);

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
%  Typical suggestions in the literature are on the form
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

% This is done within the integration-point loop.
tau=SUPGtau(CtrlVar,speed0,l,dt,CtrlVar.uvh.SUPG.tau,CtrlVar.uvh.SUPG.tauMultiplier) ;
tau0=CtrlVar.SUPG.beta0*tau;




detJw=detJ*MUA.weights(Iint);
nod=MUA.nod ;

switch lower(CtrlVar.FlowApproximation)

    case "sstream"

       
        UseMex=~Ronly && CtrlVar.UseMexFiles ;

        if ~UseMex


            [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
                uvhNodalLoopSSTREAM(detJw,nod,theta,tau0,Ronly,...
                CtrlVar,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
                Deriv,fun,...
                exx,eyy,exy,exx0,eyy0,...
                dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
                ca,sa,g,dt,...
                etaint,Eint,...
                h0barr,h1barr,...
                taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
                Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
                hint,h0int,a1int,a0int,dadhint,lambda_h) ;

        else

            if CtrlVar.UseMexFilesCPUcompare
                Tx2=Tx ; Fx2=Fx ; Ty2=Ty ; Fy2=Fy ; Th2=Th; Fh2=Fh ;
                Kxu2=Kxu ;  Kxv2=Kxv ; Kyu2= Kyu;
                Kyv2=Kyv ;  Kxh2=Kxh ; Kyh2=Kyh ;
                Khu2= Khu ; Khv2=Khv ; Khh2=Khh;
            end

            CV=1; Ronly=double(Ronly) ;

            tStart = cputime;

            [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
                uvhNodalLoopSSTREAM_mex(detJw,nod,theta,tau0,Ronly,...
                CV,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
                Deriv,fun,...
                exx,eyy,exy,exx0,eyy0,...
                dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
                ca,sa,g,dt,...
                etaint,Eint,...
                h0barr,h1barr,...
                taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
                Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
                hint,h0int,a1int,a0int,dadhint,lambda_h) ;

            tEndMex = cputime - tStart ; %fprintf("    MEX file %f \n",tEndMex)

            if CtrlVar.UseMexFilesCPUcompare
                tStart = cputime;
                [Tx2,Fx2,Ty2,Fy2,Th2,Fh2,Kxu2,Kxv2,Kyu2,Kyv2,Kxh2,Kyh2,Khu2,Khv2,Khh2]=...
                    uvhNodalLoopSSTREAM(detJw,nod,theta,tau0,Ronly,...
                    CtrlVar,Tx2,Fx2,Ty2,Fy2,Th2,Fh2,Kxu2,Kxv2,Kyu2,Kyv2,Kxh2,Kyh2,Khu2,Khv2,Khh2, ...
                    Deriv,fun,...
                    exx,eyy,exy,exx0,eyy0,...
                    dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
                    ca,sa,g,dt,...
                    etaint,Eint,...
                    h0barr,h1barr,...
                    taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
                    Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
                    hint,h0int,a1int,a0int,dadhint,lambda_h) ;

                tEndm = cputime - tStart ; fprintf(" Mex %f sec \t   m-file %f sec \t m/Mex=%f \n",tEndMex,tEndm,tEndm/tEndMex)
                fprintf("Difference between solutions: %f %f %f %f \n",norm(Tx2-Tx),norm(Fx2-Fx),norm(Ty2-Ty),norm(Fy2-Fy))
            end

        end




    case "hybrid"


        [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
            uvhNodalLoopHybrid(detJw,nod,theta,tau0,Ronly,...
            CtrlVar,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
            Deriv,fun,...
            exx,eyy,exy,exx0,eyy0,...
            dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
            ca,sa,g,dt,...
            etaint,Eint,...
            h0barr,h1barr,...
            taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
            Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
            hint,h0int,a1int,a0int,dadhint,lambda_h) ;


    case "sstream-rho"

        [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
            uvhNodalLoopSSTREAMrho(detJw,nod,theta,tau0,Ronly,...
            CtrlVar,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
            Deriv,fun,...
            exx,eyy,exy,exx0,eyy0,...
            dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
            ca,sa,g,dt,...
            etaint,Eint,...
            h0barr,h1barr,...
            taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
            Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
            hint,h0int,a1int,a0int,dadhint,lambda_h) ;





    case "sstreamTest"
        error('not finalized')

        [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
            uvhNodalLoopSSTREAMtest(detJw,nod,theta,tau0,Ronly,...
            CtrlVar,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
            Deriv,fun,...
            exx,eyy,exy,exx0,eyy0,...
            dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
            ca,sa,g,dt,...
            etaint,Eint,...
            h0barr,h1barr,...
            taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
            Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
            hint,h0int,a1int,a0int,dadhint,lambda_h) ;


    otherwise
        error("What case?")
end


end
