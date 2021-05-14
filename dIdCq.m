function dIdC=dIdCq(CtrlVar,UserVar,MUA,F,uAdjoint,vAdjoint,Meas)

narginchk(7,7)

%
% Calculates the product: dFuv/dC  \lambda
%
% If we write:
%
%   dJ/dC  = dFuv/dC lambda  + dJ/dC
%          = dI/dC + dJ/dC
%
% then this is equal to dI/dC, hence the name, although this should maybe be referred to as dIdp
%
% dI/dC is the 'misfit' contribution to dJ/dC
%
% dIdC=dIdCq(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)
%
%
%
%  < (dF/dC)^* l | phi >
%
% Example:  F = C u
%
% < u dC l | phi >
%
% if dC=phi then : <u l phi_q | phi >
%
% nodal based gradient

ndim=2;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if CtrlVar.SlidingLaw=="Budd"
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    qnod=mnod*0 ;  % just to avoid asking this again within a loop
end

if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    muknod=mnod*0 ;
end


Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hfnod=F.rhow*(Snod-Bnod)./rhonod;

uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);


% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    
    
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    qint=qnod*fun;
    mukint=muknod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;
    uAdjointint=uAdjointnod*fun;
    vAdjointint=vAdjointnod*fun;
    hfint=hfnod*fun;
    %hfint=(Sint-Bint)*F.rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    %
    % dF/dC=dtaux/dC uAdjoint + dtauy/dC vAdjoint
    %
    % dtaux/dC= He u * dbeta2/dC
    %
    % beta2= (C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    %
    
    % setting this CtrlVar field to true ensures that BasalDrag.m returns the (point) derivative
    CtrlVar.Inverse.dFuvdClambda=true;
    Ctemp= ...
        BasalDrag(CtrlVar,MUA,Heint,[],hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,[],[],[],[],[],[],[],[],qint,F.g,mukint);
    CtrlVar.Inverse.dFuvdClambda=false;
    
    if ~CtrlVar.DevelopmentVersion
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            Ctemp=log(10)*Cint.*Ctemp;
        end
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        
        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*uAdjointint+vint.*vAdjointint).*fun(Inod).*detJw;
        
    end
end

dIdCtemp=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdCtemp=dIdCtemp+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

if CtrlVar.DevelopmentVersion
    % change of variables should be done on nodal values!
    % I learned this the hard way by doing extensive tests on dJ/dgamma
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        dIdCtemp=log(10)*F.C.*dIdCtemp;
    end
end

dIdC=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,[],dIdCtemp);



end





