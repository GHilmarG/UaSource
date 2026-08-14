


function [Psid2FdAdu,Psid2FdAdv]=Psi_d2FdAdq_xi(CtrlVar,MUA,F,Psi_x,Psi_y)

narginchk(8,8)




if ~contains(CtrlVar.Inverse.InvertFor,"logAGlen",IgnoreCase=true)
    Psid2FdAdu=[];
    Psid2FdAdv=[];
    return
end


if CtrlVar.SlidingLaw~="Weertman"

    error("Psi_d2Fdpdq_xi:NotImplemented","only implemented for Weertman sliding law.")

end

if contains(CtrlVar.Inverse.Measurements,"-dhdt-")

    error("Psi_d2FdCdq_xi:NotImplemented","not implemented for dhdt meas")

end

ndim=2;
nNodes=MUA.Nnodes ;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_xnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_ynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);



Hu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Hv=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);

    h=hnod*fun;


    dudx=zeros(MUA.Nele,1);
    dudy=zeros(MUA.Nele,1);

    dvdx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);

    dPsixdx=zeros(MUA.Nele,1);
    dPsiydx=zeros(MUA.Nele,1);
    dPsixdy=zeros(MUA.Nele,1);
    dPsiydy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod) ;
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod) ;

        dPsixdx=dPsixdx+Deriv(:,1,Inod).*Psi_xnod(:,Inod);
        dPsiydx=dPsiydx+Deriv(:,1,Inod).*Psi_ynod(:,Inod);
        dPsixdy=dPsixdy+Deriv(:,2,Inod).*Psi_xnod(:,Inod);
        dPsiydy=dPsiydy+Deriv(:,2,Inod).*Psi_ynod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);
    
  
    eEff=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));
    eta0=CtrlVar.etaZero; 
    eta=real(0.5*A.^(-1./n).*eEff.^((1-n)./n))+eta0 ; 


    E1_uv = 4.*dudx+2.*dvdy;
    E2_uv = dudy+dvdx;
    F1_uv = 4.*dvdy+2.*dudx;

    E1_PsixPsiy = 4.*dPsixdx+2.*dPsiydy;
    E2_PsixPsiy = dPsixdy+dPsiydx;
    F1_PsixPsiy = 4.*dPsiydy+2.*dPsixdx;

    S = -log(10) ./ n .* (eta - eta0);
    Theta = E1_uv.*dPsixdx + E2_uv.*dPsixdy + F1_uv.*dPsiydy + E2_uv.*dPsiydx;
    R = -(1-n) .* log(10) ./ (8*n.^2) .* A.^(-1./n) .* eEff.^((1-3*n)./n);



    K1 = h .* ( R.*Theta.*E1_uv + S.*E1_PsixPsiy );
    K2 = h .* ( R.*Theta.*E2_uv + S.*E2_PsixPsiy );
    K3 = h .* ( R.*Theta.*F1_uv + S.*F1_PsixPsiy );

    for Inod=1:MUA.nod
        for Lnod=1:MUA.nod

      
      
            phi_l=fun(Lnod);
           % phi_i=fun(Inod) ;
           
           % dphldx_l=Deriv(:,1,Lnod);
           % dphidy_l=Deriv(:,2,Lnod);
            dphidx_i=Deriv(:,1,Inod);
            dphidy_i=Deriv(:,2,Inod);
            
            Hu(:,Inod,Lnod)=Hu(:,Inod,Lnod) + phi_l .* (K1.*dphidx_i + K2.*dphidy_i) .*detJw;
            Hv(:,Inod,Lnod)=Hv(:,Inod,Lnod) + phi_l .* (K2.*dphidx_i + K3.*dphidy_i) .*detJw;
            
          

        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
HuVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
HvVal=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Lnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Lnod);
        HuVal(istak+1:istak+MUA.Nele)=Hu(:,Inod,Lnod);
        HvVal(istak+1:istak+MUA.Nele)=Hv(:,Inod,Lnod);


        istak=istak+MUA.Nele;

    end
end

Psid2FdAdu=sparseUA(Iind,Jind,HuVal,nNodes,nNodes) ;  % this will be full matrix,
Psid2FdAdv=sparseUA(Iind,Jind,HuVal,nNodes,nNodes) ;  % this will be full matrix,


return

% For A/C mixing, must do this outside together with the C matrix
% 
% Psi_d2FduA_dudA=Psid2FdAdu*KdudA;  % this will be full matrix,
% Psi_d2FdvA_dvdA=Psid2FdAdv*KdvdA;
% 
% Psi_d2FduA_dudA=full(Psi_d2FduA_dudA);
% Psi_d2FdvA_dvdA=full(Psi_d2FdvA_dvdA);
% 
% Psi_d2FduA_dudA=Psi_d2FduA_dudA+Psi_d2FduA_dudA';
% Psi_d2FdvA_dvdA=Psi_d2FdvA_dvdA+Psi_d2FdvA_dvdA';
% 
% Psi_d2FdAq=Psi_d2FduA_dudA+Psi_d2FdvA_dvdA;


end







