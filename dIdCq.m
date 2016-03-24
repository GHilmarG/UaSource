function dIdC=dIdCq(CtrlVar,MUA,lx,ly,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

%function [dIdC]=dIdCq(u,v,lx,ly,S,B,h,MUA.connectivity,coordinates,nip,C,m,rho,rhow,CtrlVar)

% Here I evaluate all fields at integration points within the
% integration loop
%




ndim=2;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
lxnod=reshape(lx(MUA.connectivity,1),MUA.Nele,MUA.nod);
lynod=reshape(ly(MUA.connectivity,1),MUA.Nele,MUA.nod);


[points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    
    
    
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    rhoint=rhonod*fun;
    lxint=lxnod*fun;
    lyint=lynod*fun;
    hfint=(Sint-Bint)*rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    Ctemp= (1/m)*Heint.*(Cint+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
    
    detJw=detJ*weights(Iint);
    for Inod=1:MUA.nod
        
        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*lxint+vint.*lyint).*fun(Inod).*detJw;
        
    end
end

dIdC=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdC=dIdC+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

if CtrlVar.MeshIndependentAdjointGradients
    
    M=MassMatrix2D1dof(MUA);
    dIdCm=M\dIdC;
    
%     figure
%     subplot(1,2,1); PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; title('dIdC nonscaled')
%     subplot(1,2,2); PlotMeshScalarVariable(CtrlVar,MUA,dIdCm) ; title('M\\dIdC')
    
    dIdC=dIdCm;
    
end

end



