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

[points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    
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
    
    
    dudx=zeros(MUA.Nele,1); dvdx=zeros(MUA.Nele,1); dudy=zeros(MUA.Nele,1); dvdy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
         
        dlxdx=dlxdx+Deriv(:,1,Inod).*uAdjointnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*vAdjointnod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*uAdjointnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*vAdjointnod(:,Inod);
        
    end
    

    detJw=detJ*weights(Iint);

    exy=(dudy+dvdx)/2;    
    [~,~,~,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,dudx,dvdy,exy);
    dEtadA=dEtadA.*hint; 

    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            
            dEtadA=log(10)*AGlenInt.*dEtadA;
            
    end
    
    
    for Inod=1:MUA.nod
        T(:,Inod)=T(:,Inod)...
            -dEtadA.*((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw;
    end
end

dIdA=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdA=dIdA+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end



switch CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdA) ; title('dIdA Mesh Dependend')
        end
        
        dIACnorm=norm(dIdA);
        dIdA=MUA.M\dIdA;
        dIdA=dIdA*dIACnorm/norm(dIdA);
        
        
        if CtrlVar.Inverse.InfoLevel>=10
            fprintf('Making dIdA mesh independent by premultiplying with the inverse of the mass matrix.\n')
        end
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdA) ; title('dIdA Mesh Independend')
        end
end



end






