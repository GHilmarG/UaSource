function dIdC=dIdCq(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

% nodal based gradient

ndim=2;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(m(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);


[points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    detJ=MUA.DetJ(:,Iint);
    
    
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    rhoint=rhonod*fun;
    uAdjointint=uAdjointnod*fun;
    vAdjointint=vAdjointnod*fun;
    hfint=(Sint-Bint)*rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    %
    % dF/dC=dtaux/dC uAdjoint + dtauy/dC vAdjoint
    %
    % dtaux/dC= He u * dbeta2/dC 
    %
    % beta2= (C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    %
    Ctemp= (1./mint).*Heint.*(Cint+CtrlVar.Czero).^(-1./mint-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1./mint-1) ;
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        Ctemp=log(10)*Cint.*Ctemp;
    end
    
    detJw=detJ*weights(Iint);
    for Inod=1:MUA.nod
        
        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*uAdjointint+vint.*vAdjointint).*fun(Inod).*detJw;
        
    end
end

dIdC=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdC=dIdC+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


switch CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; title('dIdC Mesh Dependend')
        end
        
        dIdCnorm=norm(dIdC);
        dIdC=MUA.M\dIdC;
        dIdC=dIdC*dIdCnorm/norm(dIdC);
        
        
        if CtrlVar.Inverse.InfoLevel>=10
            fprintf('Making dIdC mesh independent by premultiplying with the inverse of the mass matrix.\n')
        end
        
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; title('dIdC Mesh Independend')
        end
end



end





