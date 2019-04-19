function dIdb=dIdbq(CtrlVar,MUA,uAdjoint,vAdjoint,F,dhdtres,dhdtErr)

%
% nodal-based gradients
%

[Areas,xEle,yEle,Area]=TriAreaFE(MUA.coordinates,MUA.connectivity);  % add areas a N to MUA

ndim=2;


if ~CtrlVar.IncludeMelangeModelPhysics
    uoint=[];
    voint=[];
    Coint=[];
    moint=[];
    uaint=[];
    vaint=[];
    Caint=[];
    maint=[];
else
    error('Inversion with MelangeModelPhysics not implemented.\n')
end

ca=cos(F.alpha); sa=sin(F.alpha);

if ~isempty(dhdtres)
    dhdtresnod=reshape(dhdtres(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    dhdtresnod=0;
    dhdtErrnod=0;
end

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
%bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);


rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);
    
    %[Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    
    
    hint=hnod*fun;
    %bint=bnod*fun;
    %sint=snod*fun;
    dhdtresint=dhdtresnod*fun;
    dhdtErrint=dhdtErrnod*fun;
    
    uint=unod*fun;
    vint=vnod*fun;
    
    
    rhoint=rhonod*fun;
    nint=nnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    
    
    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    uAdjointint=uAdjointnod*fun;
    vAdjointint=vAdjointnod*fun;
    Hint=Sint-Bint;
    hfint=F.rhow*Hint./rhoint;
    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    dudx=zeros(MUA.Nele,1);     dvdx=zeros(MUA.Nele,1);     dudy=zeros(MUA.Nele,1);     dvdy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1);    dlydx=zeros(MUA.Nele,1);    dlxdy=zeros(MUA.Nele,1);    dlydy=zeros(MUA.Nele,1);
    dsdx=zeros(MUA.Nele,1);     dsdy=zeros(MUA.Nele,1);
    dhdx=zeros(MUA.Nele,1);     dhdy=zeros(MUA.Nele,1);
    drhodx=zeros(MUA.Nele,1);   drhody=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
        
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
        drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);
        
        
        dlxdx=dlxdx+Deriv(:,1,Inod).*uAdjointnod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*uAdjointnod(:,Inod);
        
        dlydx=dlydx+Deriv(:,1,Inod).*vAdjointnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*vAdjointnod(:,Inod);
        
    end
    
    
    detJw=detJ*weights(Iint);
    
    exy=(dudy+dvdx)/2;
    
    
    %
    [~,~,~,~,~,~,dtaubxdh,dtaubydh] = BasalDrag(CtrlVar,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint);
    etaint=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,dudx,dvdy,exy);
    
    
    for Inod=1:MUA.nod
        
        T(:,Inod)=T(:,Inod)...
            +etaint.*((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw...
            +(dtaubxdh.*uAdjointint+dtaubydh.*vAdjointint).*fun(Inod).*detJw...
            +F.g*ca*((rhoint.*dsdx+hint.*drhodx).*uAdjointint+(rhoint.*dsdy+hint.*drhody).*vAdjointint).*fun(Inod).*detJw...
            -F.g*sa*((rhoint+hint.*drhodx).*uAdjointint).*fun(Inod).*detJw ;
        
    end
    
    if contains(CtrlVar.Inverse.Measurements,'-dhdt-')
        for Inod=1:MUA.nod
            
            dbI=dhdtresint...
                .*(dudx.*fun(Inod)+uint.*Deriv(:,1,Inod)+dvdy.*fun(Inod)+vint.*Deriv(:,2,Inod))...
                .*detJw./dhdtErrint/Area ;
            T(:,Inod)=T(:,Inod)+dbI ;
            
        end
        
    end
    
end

dIdb=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdb=dIdb+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end



switch CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdb) ; title('dIdb Mesh Dependend')
        end
        
        dIdbnorm=norm(dIdb);
        dIdb=MUA.M\dIdb;
        dIdb=dIdb*dIdbnorm/norm(dIdb);
        
        
        if CtrlVar.Inverse.InfoLevel>=10
            fprintf('Making dIdb mesh independent by premultiplying with the inverse of the mass matrix.\n')
        end
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdb) ; title('dIdb Mesh Independend')
        end
end



end






