



function KFqq=Fqq(CtrlVar,MUA,F,Psi_x,Psi_y)

%%
%
% $$\langle \Psi_x \mid \delta^2_{qq}\mathcal{F}_x[\xi_{,l},\xi_{,m}] \rangle + \langle \Psi_y \mid \delta^2_{qq}\mathcal{F}_y[\xi_{,l},\xi_{,m}] \rangle $$
%
%
%%


narginchk(5,5)

ndim=2;



% eta0=CtrlVar.etaZero;
Eps0=CtrlVar.EpsZero;
C0=CtrlVar.Czero;
u0=CtrlVar.SpeedZero;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
% snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
% bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod


Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_xnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_ynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

Fuu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Fvu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Fuv=zeros(MUA.Nele,MUA.nod,MUA.nod);
Fvv=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);

    h=hnod*fun;

    % s=snod*fun;

    B=Bnod*fun;
    S=Snod*fun;
    H=S-B;

    u=unod*fun;
    v=vnod*fun;

    Psi_x=Psi_xnod*fun;
    Psi_y=Psi_ynod*fun;

    n=nnod*fun;
    A=AGlennod*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;

    m=mnod*fun;
    C=Cnod*fun;
    C(C<CtrlVar.Cmin)=CtrlVar.Cmin;

    rho=rhonod*fun;


    hf=F.rhow*H./rho;                                   % this is linear, so fine to evaluate at int in this manner
    He = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
   



    dudx=zeros(MUA.Nele,1);
    dudy=zeros(MUA.Nele,1);

    dvdx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);

    dPsi_xdx=zeros(MUA.Nele,1);
    dPsi_ydx=zeros(MUA.Nele,1);
    dPsi_xdy=zeros(MUA.Nele,1);
    dPsi_ydy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod) ;
        dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod) ;

        dPsi_xdx=dPsi_xdx+Deriv(:,1,Inod).*Psi_xnod(:,Inod);
        dPsi_ydx=dPsi_ydx+Deriv(:,1,Inod).*Psi_ynod(:,Inod);
        dPsi_xdy=dPsi_xdy+Deriv(:,2,Inod).*Psi_xnod(:,Inod);
        dPsi_ydy=dPsi_ydy+Deriv(:,2,Inod).*Psi_ynod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);

    E_1 = 4.*dudx+2.*dvdy;
    E_2 = dudy+dvdx;
    F_1 = 4.*dvdy+2.*dudx;
    e_eff = sqrt((dudx).^2+(dvdy).^2+dudx .* dvdy+0.25.*(dvdx+dudy).^2+Eps0.^2);
    u_e = sqrt(u.^2+v.^2+u0.^2);
    G=He; 


    for Inod=1:MUA.nod
        for Knod=1:MUA.nod

            phi_k=fun(Knod);
            phi_i=fun(Inod) ;
            dphidx_i=Deriv(:,1,Inod);
            dphidy_i=Deriv(:,2,Inod);
            dphidx_k=Deriv(:,1,Knod);
            dphidy_k=Deriv(:,2,Knod);

            du_eta_i = (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./(n)).* (E_1 .* dphidx_i + E_2 .* dphidy_i);
            dv_eta_i = (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./(n)).* (F_1 .* dphidy_i + E_2 .* dphidx_i);

            du_eta_k = (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./(n)).* (E_1 .* dphidx_k + E_2 .* dphidy_k);
            dv_eta_k = (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./(n)).* (F_1 .* dphidy_k + E_2 .* dphidx_k);

            du_beta2_i = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-3.*m)./(m)) .* u .*  phi_i;
            dv_beta2_i = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-3.*m)./(m)) .* v .*  phi_i;

            du_beta2_k = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-3.*m)./(m)) .* u .*  phi_k;
            dv_beta2_k = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-3.*m)./(m)) .* v .*  phi_k;

            d2uu_beta2_ik = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-5.*m)./(m)).*((1-2.*m)./(m) .* u.^2 + (v.^2+u0.^2)).*phi_i .* phi_k;
            d2vv_beta2_ik = (1-m)./(m) .* (C+C0).^(-1./m) .* u_e.^((1-5.*m)./(m)).*((1-2.*m)./(m) .* v.^2 + (u.^2+u0.^2)).*phi_i .* phi_k;
            d2uv_beta2_ik = ((1-m).*(1-3.*m))./(m.^2) .* (C+C0).^(-1./m) .* u_e.^((1-5.*m)./(m)) .* u.*v .*  phi_i .* phi_k;
            d2vu_beta2_ik = d2uv_beta2_ik;


            d2uu_eta_ik = ((1-n).*(1-3.*n))./(32.*n.^2) .* A.^(-1./n) .* e_eff.^((1-5.*n)./(n)).*(E_1 .* dphidx_i+E_2 .* dphidy_i ) .* (E_1 .* dphidx_k+E_2 .* dphidy_k) ...
                        + (1-n)./(2.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./n).*( dphidx_i .* dphidx_k+0.25.* dphidy_i .* dphidy_k);

            d2vv_eta_ik = ((1-n).*(1-3.*n))./(32.*n.^2) .* A.^(-1./n) .* e_eff.^((1-5.*n)./n).*(F_1 .* dphidy_i+E_2 .* dphidx_i) .* (F_1 .* dphidy_k+E_2 .* dphidx_k) ...
            + (1-n)./(2.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./n).*( dphidy_i .* dphidy_k+0.25.* dphidx_i .* dphidx_k) ; 


            % d2vu_eta_ik=d2uv_eta_ki i.e. we need to swap both u <-> v, and i <-> j , otherwise they are not equal!
     
            d2uv_eta_ik = ((1-n).*(1-3.*n))./(32.*n.^2) .* A.^(-1./n) .* e_eff.^((1-5.*n)./n).*(E_1 .* dphidx_i+E_2 .* dphidy_i) .* (F_1 .* dphidy_k+E_2 .* dphidx_k) ...
            + (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./n) .* (2 .* dphidx_i .* dphidy_k+ dphidy_i .* dphidx_k) ; 


            d2vu_eta_ik = ((1-n).*(1-3.*n))./(32.*n.^2) .* A.^(-1./n) .* e_eff.^((1-5.*n)./n).*(E_1 .* dphidx_k+E_2 .* dphidy_k) .* (F_1 .* dphidy_i+E_2 .* dphidx_i) ...
            + (1-n)./(8.*n) .* A.^(-1./n) .* e_eff.^((1-3.*n)./n).*(2 .* dphidx_k .* dphidy_i+ dphidy_k .* dphidx_i) ; 


            F_uu_ik = ...
                h .* d2uu_eta_ik .* (4 .* dudx+2 .* dvdy) .*   dPsi_xdx ...
            +   h .* du_eta_i .* 4 .* dphidx_k .*   dPsi_xdx ...
            +   h .* du_eta_k .* 4 .* dphidx_i .*   dPsi_xdx ...
            +   h .* d2uu_eta_ik .* ( dudy+ dvdx) .*   dPsi_xdy ...
            +   h .* du_eta_i .*  dphidy_k .*   dPsi_xdy ...
            +   h .* du_eta_k .*  dphidy_i .*   dPsi_xdy ...
            +   G .* d2uu_beta2_ik .* u .*  Psi_x ...
            +   G .* du_beta2_i .*  phi_k .*  Psi_x ...
            +   G .* du_beta2_k .*  phi_i .*  Psi_x ... 
            +   h .* d2uu_eta_ik .* (4 .* dvdy+2 .* dudx) .*   dPsi_ydy ...
            +   h .* du_eta_i .* 2 .* dphidx_k .*   dPsi_ydy ...
            +   h .* du_eta_k .* 2 .* dphidx_i .*   dPsi_ydy ...
            +   h .* d2uu_eta_ik .* ( dvdx+ dudy) .*   dPsi_ydx ...
            +   h .* du_eta_i .*  dphidy_k .*   dPsi_ydx ... 
            +   h .* du_eta_k .*  dphidy_i .*   dPsi_ydx ...
            +   G .* d2uu_beta2_ik .* v .*  Psi_y;



            F_vv_ik = ...
                h .* d2vv_eta_ik .* (4 .* dudx+2 .* dvdy) .*   dPsi_xdx ...
            +   h .* dv_eta_i .* 2 .* dphidy_k .*   dPsi_xdx ...
            +   h .* dv_eta_k .* 2 .* dphidy_i .*   dPsi_xdx ...
            +   h .* d2vv_eta_ik .* ( dudy+ dvdx) .*   dPsi_xdy ...
            +   h .* dv_eta_i .*  dphidx_k .*   dPsi_xdy ...
            +   h .* dv_eta_k .*  dphidx_i .*   dPsi_xdy ...
            +   G .* d2vv_beta2_ik .* u .*  Psi_x ...
            +   h .* d2vv_eta_ik .* (4.* dvdy+2 .* dudx) .*   dPsi_ydy ...
            +   h .* dv_eta_i .* 4 .* dphidy_k .*   dPsi_ydy ...
            +   h .* dv_eta_k .* 4 .* dphidy_i .*   dPsi_ydy ...
            +   h .* d2vv_eta_ik .* ( dvdx+ dudy) .*   dPsi_ydx ...
            +   h .* dv_eta_i .*  dphidx_k .*   dPsi_ydx ...
            +   h .* dv_eta_k .*  dphidx_i .*   dPsi_ydx ...
            +   G .* d2vv_beta2_ik .* v .*  Psi_y ...
            +   G .* dv_beta2_i .*  phi_k .*  Psi_y ... 
            +   G .* dv_beta2_k .*  phi_i .*  Psi_y ;



            F_uv_ik = ...
                h .* d2uv_eta_ik .* (4.* dudx+2 .* dvdy) .*   dPsi_xdx ...
            +   h .* du_eta_i .* 2 .* dphidy_k .*   dPsi_xdx ...
            +   h .* dv_eta_k .* 4 .* dphidx_i .*   dPsi_xdx ...
            +   h .* d2uv_eta_ik .* ( dudy+ dvdx) .*   dPsi_xdy ...
            +   h .* du_eta_i .*  dphidx_k .*   dPsi_xdy ...
            +   h .* dv_eta_k .*  dphidy_i .*   dPsi_xdy ...
            +   G .* d2uv_beta2_ik .* u .*  Psi_x ...
            +   G .* dv_beta2_k .*  phi_i .*  Psi_x ...
            +   h .* d2uv_eta_ik .* (4.* dvdy+2 .* dudx) .*   dPsi_ydy ...
            +   h .* du_eta_i .* 4 .* dphidy_k .*   dPsi_ydy ...
            +   h .* dv_eta_k .* 2 .* dphidx_i .*   dPsi_ydy ...
            +   h .* d2uv_eta_ik .* ( dvdx+ dudy) .*   dPsi_ydx ...
            +   h .* du_eta_i .*  dphidx_k .*   dPsi_ydx ...
            +   h .* dv_eta_k .*  dphidy_i .*   dPsi_ydx ...
            +   G .* d2uv_beta2_ik .* v .*  Psi_y ...
            +   G .* du_beta2_i .*  phi_k .*  Psi_y ;



            F_vu_ik = ...
                h .* d2vu_eta_ik .* (4.* dudx+2 .* dvdy) .*   dPsi_xdx ...
            +   h .* dv_eta_i .* 4 .* dphidx_k .*   dPsi_xdx ...
            +   h .* du_eta_k .* 2 .* dphidy_i .*   dPsi_xdx ...
            +   h .* d2vu_eta_ik .* ( dudy+ dvdx) .*   dPsi_xdy ...
            +   h .* dv_eta_i .*  dphidy_k .*   dPsi_xdy ...
            +   h .* du_eta_k .*  dphidx_i .*   dPsi_xdy ...
            +   G .* d2vu_beta2_ik .* u .*  Psi_x ...
            +   G .* dv_beta2_i .*  phi_k .*  Psi_x ...
            +   h .* d2vu_eta_ik .* (4 .* dvdy+2 .* dudx) .*   dPsi_ydy ...
            +   h .* dv_eta_i .* 2 .* dphidx_k .*   dPsi_ydy ...
            +   h .* du_eta_k .* 4 .* dphidy_i .*   dPsi_ydy ...
            +   h .* d2vu_eta_ik .* ( dvdx+ dudy) .*   dPsi_ydx ...
            +   h .* dv_eta_i .*  dphidy_k .*   dPsi_ydx ...
            +   h .* du_eta_k .*  dphidx_i .*   dPsi_ydx ...
            +   G .* d2vu_beta2_ik .* v .*  Psi_y ...
            +   G .* du_beta2_k .*  phi_i .*  Psi_y ; 




            Fuu(:,Inod,Knod)=Fuu(:,Inod,Knod) + F_uu_ik.*detJw;
            Fuv(:,Inod,Knod)=Fuv(:,Inod,Knod) + F_uv_ik.*detJw;
            Fvu(:,Inod,Knod)=Fvu(:,Inod,Knod) + F_vu_ik.*detJw;
            Fvv(:,Inod,Knod)=Fvv(:,Inod,Knod) + F_vv_ik.*detJw;


        end
    end
end


Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

%
%
% KFqq= [ Fuu Fuv ]
%       [ Fvu Fvv ] 
%
%
%
%

Nnodes=MUA.Nnodes;
for Inod=1:MUA.nod
    for Knod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod); 
        Xval(istak+1:istak+MUA.Nele)=Fuu(:,Inod,Knod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+Nnodes; 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod)+Nnodes; 
        Xval(istak+1:istak+MUA.Nele)=Fvv(:,Inod,Knod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod)+Nnodes; 
        Xval(istak+1:istak+MUA.Nele)=Fuv(:,Inod,Knod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+Nnodes; 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod); 
        Xval(istak+1:istak+MUA.Nele)=Fvu(:,Inod,Knod);
        istak=istak+MUA.Nele;


    end
end



KFqq=sparseUA(Iind,Jind,Xval,2*Nnodes,2*Nnodes);

KFqq=(KFqq'+KFqq)/2;





end







