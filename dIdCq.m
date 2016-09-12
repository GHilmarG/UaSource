function dIdC=dIdCq(CtrlVar,MUA,lx,ly,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

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
    mint=mnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    rhoint=rhonod*fun;
    lxint=lxnod*fun;
    lyint=lynod*fun;
    hfint=(Sint-Bint)*rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    Ctemp= (1./mint).*Heint.*(Cint+CtrlVar.CAdjointZero).^(-1./mint-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1./mint-1) ;
    
    detJw=detJ*weights(Iint);
    for Inod=1:MUA.nod
        
        T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*lxint+vint.*lyint).*fun(Inod).*detJw;
        
    end
end

dIdC=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdC=dIdC+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

if ~strcmpi(CtrlVar.MeshIndependentAdjointGradients,'I')
    
    switch upper(CtrlVar.MeshIndependentAdjointGradients)
        
        case 'P'
            P=NodalFormFunctionInfluence(MUA);
            dIdCm=dIdC./P;
        case 'M'
            M=MassMatrix2D1dof(MUA);
            dIdCm=M\dIdC;
    end
    
    if  CtrlVar.doplots && CtrlVar.InfoLevelAdjoint>100
        figure ;
        PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ;
        title('Euclidian dIdC i.e. \deltaJ(C,N) ')
        figure
        PlotMeshScalarVariable(CtrlVar,MUA,dIdCm) ;
        title(['Ritz representation (',CtrlVar.MeshIndependentAdjointGradients,')'])
    end
    
    dd=dIdCm'*dIdC/(norm(dIdCm)*norm(dIdC));
    ddAngle=acosd(dd);
    fprintf(' Angle between Euclidian and Ritz gradients is %g degrees.\n',ddAngle)
    
    dIdC=dIdCm;
    
end

%dIdC=median(C)*dIdC;  % reasonable scaling I think

end





