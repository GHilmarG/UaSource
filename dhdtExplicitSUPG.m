function [UserVar,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F,BCs)
%%
% Calculates dh/dt from flux divergence as
%
%   dh/dt = a -  ( dqx/dx + dqy/dy)
%
%   [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F)
%
% Note: 
%   -Uses u=F.ub, and hence only correct for plug flow, e.g. SSA
%   -Does not account for spatially variable density
%
% Projects the values directly onto nodes.
%
% Uses SUPG, but very doubtfull that this SUPG treatment is required. In fact I'm not really solving for dhdt, 
% I'm just evaluating the  spatial derivatives terms on the integration points, and then projecting onto the nodes.
%
% see also : ProjectFintOntoNodes
%
%

narginchk(5,5)

if nargin<5
    BCs=[];
end

% I tested various tau options agains the typical Galerkin without upwinding
% for an glacier peak example, and found taus created oscillations for 
% small velocities, i.e. too much negative diffusion.

CtrlVar.Tracer.SUPG.tau="tau2";
% CtrlVar.Tracer.SUPG.tau="taus";      % too much negative diffusion because of the tau=l/(2u)  goes to infinity as u to zero 
% CtrlVar.Tracer.SUPG.tau="tau1";    % results identical to noSUPG 
% CtrlVar.Tracer.SUPG.tau="taut";    % results identical to noSUPG 
% CtrlVar.Tracer.SUPG.tau="tau2";    % results identical to noSUPG 



ndim=2; dof=1; neq=dof*MUA.Nnodes;
EleAreas=[];
tauSUPG=CalcSUPGtau(CtrlVar,EleAreas,F.ub,F.vb,CtrlVar.dt,MUA);


anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

tauSUPGnod=reshape(tauSUPG(MUA.connectivity,1),MUA.Nele,MUA.nod);

% [points,weights]=sample('triangle',MUA.nip,ndim);

dd=zeros(MUA.Nele,MUA.nod,MUA.nod);
b=zeros(MUA.Nele,MUA.nod);


% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    aint=anod*fun;
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    
    tauSUPGint=tauSUPGnod*fun;
    

    
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    
    for Inod=1:MUA.nod
        
        SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod));
        SUPGdetJw=SUPG.*detJw;
        
        for Jnod=1:MUA.nod
            
            dd(:,Inod,Jnod)=dd(:,Inod,Jnod)+fun(Jnod).*SUPGdetJw;
            
        end
        
        
        tx=(dhdx.*uint+hint.*dudx);
        ty=(dhdy.*vint+hint.*dvdy);
        
        %term=(aint-tx-ty).*fun(Inod).*detJw;
        term=(aint-tx-ty).*SUPGdetJw;
        
        b(:,Inod)=b(:,Inod)+term;
        
        
    end
end

% assemble right-hand side

rh=sparseUA(neq,1);
for Inod=1:MUA.nod
    rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b(:,Inod),neq,1);
end

if ~isfield(MUA,'M')
    MUA.M=MassMatrix2D1dof(MUA);
end


Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=dd(:,Inod,Jnod);
        istak=istak+MUA.Nele;
    end
end

Msupg=sparseUA(Iind,Jind,Xval,neq,neq);
[hL,hRhs]=createLh(MUA.Nnodes,BCs.dhdtFixedNode,BCs.dhdtFixedValue,BCs.dhdtTiedNodeA,BCs.dhdtTiedNodeB);

% CtrlVar.SymmSolver='AugmentedLagrangian';
x0=zeros(MUA.Nnodes,1) ; y0=hRhs*0;


[dhdt,dhdtlambda]=solveKApe(Msupg,hL,rh,hRhs,x0,y0,CtrlVar);
dhdt=full(dhdt);


% Now there is an issue here regarding what to do about dhdt<0 when h<=thickmin

I=(F.h<=CtrlVar.ThickMin)  & (dhdt< 0) ;
dhdt(I)=0 ; 



end