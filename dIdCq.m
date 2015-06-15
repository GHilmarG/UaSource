function [dIdC]=dIdCq(u,v,lx,ly,S,B,h,connectivity,coordinates,nip,C,m,rho,rhow,CtrlVar)
    
    
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2;
    
    
    hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
    unod=reshape(u(connectivity,1),Nele,nod);
    vnod=reshape(v(connectivity,1),Nele,nod);
    Cnod=reshape(C(connectivity,1),Nele,nod);
    Bnod=reshape(B(connectivity,1),Nele,nod);
    Snod=reshape(S(connectivity,1),Nele,nod);
    lxnod=reshape(lx(connectivity,1),Nele,nod);
    lynod=reshape(ly(connectivity,1),Nele,nod);
    
    
    [points,weights]=sample('triangle',nip,ndim);
    
    
    T=zeros(Nele,nod);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        hint=hnod*fun;
        uint=unod*fun;
        vint=vnod*fun;
        Cint=Cnod*fun;
        Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible that Cint is less than any of the nodal values
        Bint=Bnod*fun;
        Sint=Snod*fun;
        lxint=lxnod*fun;
        lyint=lynod*fun;
        
        hfint=(Sint-Bint)*rhow/rho;
        kH=CtrlVar.kH;
        Heint = HeavisideApprox(kH,hint-hfint);
        detJw=detJ*weights(Iint);
        
        Ctemp= (1/m)*Heint.*(Cint+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
        for Inod=1:nod
            
            T(:,Inod)=T(:,Inod)+Ctemp.*(uint.*lxint+vint.*lyint).*fun(Inod).*detJw;
            
        end
    end
    
    dIdC=zeros(Nnodes,1);
    
    for Inod=1:nod
        dIdC=dIdC+sparse(connectivity(:,Inod),ones(Nele,1),T(:,Inod),Nnodes,1);
    end
    
end



