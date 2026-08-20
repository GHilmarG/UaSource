


function [KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%%
%
% Calculates
%
%
% $$
% \mathcal{F}^{pq}_{Au,li} =  \langle \Psi_x , \delta^2_{Au}\mathcal{F}_x[\phi_l,\phi_i] \rangle   + \langle \Psi_y,\delta^2_{Au}\mathcal{F}_y[\phi_l,\phi_i] \rangle
% $$
%
% $$ 
% \mathcal{F}^{pq}_{Av,li} = \langle \Psi_x , \delta^2_{Av}\mathcal{F}_x[\phi_l,\phi_i] \rangle  + \langle \Psi_y,\delta^2_{Av}\mathcal{F}_y[\phi_l,\phi_i] \rangle
% $$
%
%


if ~contains(CtrlVar.Inverse.InvertFor,"logAGlen",IgnoreCase=true)
    KFAu=[];
    KFAv=[];
    return
end


if CtrlVar.SlidingLaw~="Weertman"

    error("FAuv:NotImplemented","only implemented for Weertman sliding law.")

end

if contains(CtrlVar.Inverse.Measurements,"-dhdt-")

    error("FAuv:NotImplemented","not implemented for dhdt meas")

end

ndim=2;
nNodes=MUA.Nnodes ;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_xnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_ynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);



Hu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Hv=zeros(MUA.Nele,MUA.nod,MUA.nod);

Eps0=CtrlVar.EpsZero; 
eta0=CtrlVar.etaZero ; 

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);

    h=hnod*fun;
    n=nnod*fun;
    A=AGlennod*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;

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
    
  
    %eEff=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));
   
    eEff = sqrt((dudx).^2+(dvdy).^2+dudx .* dvdy+0.25.*(dvdx+dudy).^2+Eps0.^2);
    
    eta=real(0.5*A.^(-1./n).*eEff.^((1-n)./n))+eta0 ; 


    E1_uv = 4.*dudx+2.*dvdy;
    E2_uv = dudy+dvdx;
    F1_uv = 4.*dvdy+2.*dudx;

    E1_PsixPsiy = 4.*dPsixdx+2.*dPsiydy;
    E2_PsixPsiy = dPsixdy+dPsiydx;
    F1_PsixPsiy = 4.*dPsiydy+2.*dPsixdx;

    % RAW (A-space, no log/A-scaling) building blocks -- Stilde is exactly EffectiveViscositySSTREAM's detadA
    Stilde = -A.^(-1./n-1) .* eEff.^((1-n)./n) ./ (2*n);
    Rtilde = -(1-n) ./ (8*n.^2) .* A.^(-1./n-1) .* eEff.^((1-3*n)./n);

    Theta = E1_uv.*dPsixdx + E2_uv.*dPsixdy + F1_uv.*dPsiydy + E2_uv.*dPsiydx;

    K1 = h .* ( Rtilde.*Theta.*E1_uv + Stilde.*E1_PsixPsiy );
    K2 = h .* ( Rtilde.*Theta.*E2_uv + Stilde.*E2_PsixPsiy );
    K3 = h .* ( Rtilde.*Theta.*F1_uv + Stilde.*F1_PsixPsiy );


    % 
    % 
    % 
    % S = -log(10) ./ n .* (eta - eta0);
    % Theta = E1_uv.*dPsixdx + E2_uv.*dPsixdy + F1_uv.*dPsiydy + E2_uv.*dPsiydx;
    % R = -(1-n) .* log(10) ./ (8*n.^2) .* A.^(-1./n) .* eEff.^((1-3*n)./n);
    % 
    % 
    % 
    % K1 = h .* ( R.*Theta.*E1_uv + S.*E1_PsixPsiy );
    % K2 = h .* ( R.*Theta.*E2_uv + S.*E2_PsixPsiy );
    % K3 = h .* ( R.*Theta.*F1_uv + S.*F1_PsixPsiy );

    for Inod=1:MUA.nod
        for Lnod=1:MUA.nod

            
            phi_i=fun(Inod);

            dphidx_l=Deriv(:,1,Lnod);
            dphidy_l=Deriv(:,2,Lnod);

            Hu(:,Inod,Lnod)=Hu(:,Inod,Lnod) - phi_i .* (K1.*dphidx_l + K2.*dphidy_l) .*detJw;
            Hv(:,Inod,Lnod)=Hv(:,Inod,Lnod) - phi_i .* (K2.*dphidx_l + K3.*dphidy_l) .*detJw;
            
          

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

KFAu=sparse(Iind,Jind,HuVal,nNodes,nNodes) ;  
KFAv=sparse(Iind,Jind,HvVal,nNodes,nNodes) ;  

ln10 = log(10);
D_A = spdiags(F.AGlen(:)*ln10, 0, nNodes, nNodes);   % nodal A values, NOT AGlennod

KFAu = D_A * KFAu;
KFAv = D_A * KFAv;





if CtrlVar.Inverse.TestDirectAdjoint.isTrue


    iColumn=randi(MUA.Nnodes);  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.

    CtrlVar.log10Derivatives=true;


    %% FCu test
    u0 = F.ub;  
    v0 = F.vb;
    hstep = 1e-6*max(abs(u0)+abs(v0));

    F.ub = u0; F.ub(iColumn) = F.ub(iColumn) - hstep;
    bMinus = dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.ub = u0; F.ub(iColumn) = F.ub(iColumn) + hstep;
    bPlus = dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.ub = u0;

    HAu_col_FD = (bPlus - bMinus)/(2*hstep);
    Diff=norm(KFAu(:,iColumn) - HAu_col_FD)/norm(HAu_col_FD);
    fprintf("FAu: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)


    % FCv test
  
  

    F.vb = v0; F.vb(iColumn) = F.vb(iColumn) - hstep;
    bMinus = dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.vb = v0; F.vb(iColumn) = F.vb(iColumn) + hstep;
    bPlus = dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.vb = v0;

    HAv_col_FD = (bPlus - bMinus)/(2*hstep);
    Diff=norm(KFAv(:,iColumn) - HAv_col_FD)/norm(HAv_col_FD);
    fprintf("FAv: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

    FAuTest=FindOrCreateFigure("Test: FAu") ; 

    plot(KFAu(:,iColumn),HAu_col_FD,"or") ; axis equal ;
    hold on ;
    plot([min(HAu_col_FD) max(HAu_col_FD)],[min(HAu_col_FD) max(HAu_col_FD)],"--k")

    % ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

    xlabel("Direct-Adjoint",Interpreter="latex")  ;
    ylabel("Finite difference",Interpreter="latex")
    title("$\langle \Psi_x | \delta^2_{v_x A} F_x \rangle $",Interpreter="latex")
    subtitle("Comparison is here for one random column",Interpreter="latex")


    FAvTest=FindOrCreateFigure("Test: FAv") ; 

    plot(KFAv(:,iColumn),HAv_col_FD,"or") ; axis equal ;
    hold on ;
    plot([min(HAv_col_FD) max(HAv_col_FD)],[min(HAv_col_FD) max(HAv_col_FD)],"--k")

    % ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

    xlabel("Direct-Adjoint",Interpreter="latex")  ;
    ylabel("Finite difference",Interpreter="latex")
    title("$\langle  \Psi_y | \delta^2_{v_y A} F_y   \rangle $",Interpreter="latex")
    subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")
    drawnow
 


    fprintf("FAuv: Inspect in debugger and then continue: [F5]\n")
    keyboard

end

end







