function dIdA=dIdAq(CtrlVar,MUA,lx,ly,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

%
% nodal-based gradients 
%

ndim=2;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
lxnod=reshape(lx(MUA.connectivity,1),MUA.Nele,MUA.nod);
lynod=reshape(ly(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);

[~,~,~,~,~,~,~,e]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; 
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        detJ=MUA.DetJ(:,Iint);
        Deriv=MUA.Deriv(:,:,:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    hint=hnod*fun;
    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    
    
    dudx=zeros(MUA.Nele,1); dvdx=zeros(MUA.Nele,1); dudy=zeros(MUA.Nele,1); dvdy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
        dlxdx=dlxdx+Deriv(:,1,Inod).*lxnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*lynod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*lxnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*lynod(:,Inod);
        
    end
    
    detJw=detJ*weights(Iint);
    %dEtadA=-real(hint.*AGlenInt.^(-1/n-1).*e(:,Iint).^((1-n)/n))/(2*n);
    dEtadA=-real(hint.*(AGlenInt+CtrlVar.AGlenAdjointZero).^(-1/n-1).*(e(:,Iint)+CtrlVar.AdjointEpsZero).^((1-n)/n))/(2*n);
    
    for Inod=1:MUA.nod
        T(:,Inod)=T(:,Inod)...
            -dEtadA.*((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw;
    end
end

dIdA=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdA=dIdA+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


if ~strcmpi(CtrlVar.MeshIndependentAdjointGradients,'I')
    
    switch upper(CtrlVar.MeshIndependentAdjointGradients)
        
        case 'P'
            P=NodalFormFunctionInfluence(MUA);
            dIdAm=dIdA./P;
        case 'M'
            M=MassMatrix2D1dof(MUA);
            dIdAm=M\dIdA;
    end
    
    if  CtrlVar.doplots && CtrlVar.InfoLevelAdjoint>100
        figure ;
        PlotMeshScalarVariable(CtrlVar,MUA,dIdA) ; 
        title('dIdA nonscaled, i.e. \deltaJ(A,N) ')
        figure
        PlotMeshScalarVariable(CtrlVar,MUA,dIdAm) ; 
        title(['Mesh-independent representation (',CtrlVar.MeshIndependentAdjointGradients,')'])
    end
    
    dIdA=dIdAm;
    
end


end



