function dIdA=dIdAq(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

%
% nodal-based gradients
%

ndim=2;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);

% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        detJ=MUA.DetJ(:,Iint);
        Deriv=MUA.Deriv(:,:,:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    hint=hnod*fun;
    nint=nnod*fun;
    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    
    
    exx=zeros(MUA.Nele,1); exy=zeros(MUA.Nele,1);  eyy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod)); 
        
        
        dlxdx=dlxdx+Deriv(:,1,Inod).*uAdjointnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*vAdjointnod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*uAdjointnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*vAdjointnod(:,Inod);
        
    end
    

    detJw=detJ*MUA.weights(Iint);

    
    [~,~,~,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);
    %dEtadA=dEtadA.*hint;
    
    if ~CtrlVar.DevelopmentVersion
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            
            dEtadA=log(10)*AGlenInt.*dEtadA;
            
        end
    end
    
    for Inod=1:MUA.nod
        T(:,Inod)=T(:,Inod)...
            -dEtadA.*hint.*((4*exx+2*eyy).*dlxdx+2*exy.*dlxdy+(4*eyy+2*exx).*dlydy+2*exy.*dlydx).*fun(Inod).*detJw;
        %-dEtadA.*((4*exx+2*eyy).*dlxdx+(dudy+dvdx).*dlxdy+(4*eyy+2*exx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw;
        
 
        
    end
end

dIdA=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdA=dIdA+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

if CtrlVar.DevelopmentVersion
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        
        dIdA=log(10)*AGlen.*dIdA;
        
    end
end


Happrox=MUA.M/MUA.Area;
dIdA=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,Happrox,dIdA);

end






