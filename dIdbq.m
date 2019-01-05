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
bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
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
    bint=bnod*fun;
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
    Hint=Sint-Bint;
    
    
    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    uAdjointint=uAdjointnod*fun;
    vAdjointint=vAdjointnod*fun;
    
    hfint=F.rhow*Hint./rhoint;
    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    HeHint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0);
    deltaHint=DiracDelta(CtrlVar.kH,Hint,CtrlVar.Hh0);
    dint = HeHint.*(Sint-bint);  % draft
    
    
    dudx=zeros(MUA.Nele,1);     dvdx=zeros(MUA.Nele,1);     dudy=zeros(MUA.Nele,1);     dvdy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1);    dlydx=zeros(MUA.Nele,1);    dlxdy=zeros(MUA.Nele,1);    dlydy=zeros(MUA.Nele,1);
    dsdx=zeros(MUA.Nele,1);     dsdy=zeros(MUA.Nele,1);
    dbdx=zeros(MUA.Nele,1);     dbdy=zeros(MUA.Nele,1);
    dhdx=zeros(MUA.Nele,1);     dhdy=zeros(MUA.Nele,1);
    drhodx=zeros(MUA.Nele,1);   drhody=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
        
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        
        dbdx=dbdx+Deriv(:,1,Inod).*bnod(:,Inod);
        dbdy=dbdy+Deriv(:,2,Inod).*bnod(:,Inod);
        
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
    exx=dudx;
    eyy=dvdy;
    
    %
    [~,~,~,~,~,~,dtaubxdh,dtaubydh] = BasalDrag(CtrlVar,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint);
    etaint=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,dudx,dvdy,exy);
    
    
    for Inod=1:MUA.nod
        
        t1=-ca*F.g*((rhoint.*hint-F.rhow*dint).*Deriv(:,1,Inod)-(rhoint.*fun(Inod)-F.rhow*HeHint.*fun(Inod)).*dbdx).*uAdjointint ...
            -rhoint.*F.g.*fun(Inod).*sa.*uAdjointint;
        t2=ca*F.g.*(-rhoint.*hint.*fun(Inod)+F.rhow.*dint.*(HeHint+deltaHint.*(Sint-bint)).*fun(Inod)).*dlxdx;
        
        t3=-fun(Inod).*etaint.*(4*exx+2*eyy).*dlxdx;
        t4=-fun(Inod).*etaint.*2.*exy.*dlxdy;
        t5=-dtaubxdh.*uAdjointint.*fun(Inod);
       
        
        Fx=(t1+t2).*detJw;
        Tx=(t3+t4+t5).*detJw;
        
        
        t1=-ca*F.g* ((rhoint.*hint-F.rhow*dint).*Deriv(:,2,Inod)-(rhoint.*fun(Inod)-F.rhow*HeHint.*fun(Inod)).*dbdy).*vAdjointint; % t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod)*ca;
        t2=ca*F.g.*(-rhoint.*hint.*fun(Inod)+F.rhow.*dint.*(HeHint+deltaHint.*(Sint-bint)).*fun(Inod)  ).*dlydy ; % t2=0.5*ca*g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
        
        t3=-fun(Inod).*etaint.*(4*eyy+2*exx).*dlydy; % t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
        t4=-fun(Inod).*etaint.*2.*exy.*dlydx ; % t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
        t5=-dtaubydh.*vAdjointint.*fun(Inod);   % 5=tauy.*fun(Inod);
        
        
        
        Fy=(t1+t2).*detJw;
        Ty=(t3+t4+t5).*detJw;
        
        
        T(:,Inod)=T(:,Inod)-Tx+Fx-Ty+Fy;   % opposite sign to K because of the Newton sign
        
    end
    
    
    
    % If dh/dt=d(s-b)/dt are measurments, and the control parameter p=b
    % then dJ/db=dI/db + dR/db,
    % has contributions from dI/db and dR/db
    % rather than just the R term otherwise often is the
    % case for control variables.
    %
    % calcualte: dI/db
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






