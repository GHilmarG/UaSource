function KFqq=Fqq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

%%
%
% $$\langle \Psi_x \mid \delta^2_{qq}\mathcal{F}_x[\xi_{,l},\xi_{,m}] \rangle + \langle \Psi_y \mid \delta^2_{qq}\mathcal{F}_y[\xi_{,l},\xi_{,m}] \rangle $$
%
%%

narginchk(7,7)

ndim=2;

Eps0=CtrlVar.EpsZero;
C0=CtrlVar.Czero;
u0=CtrlVar.SpeedZero;

Nele=MUA.Nele;
nod=MUA.nod;

hnod=reshape(F.h(MUA.connectivity,1),Nele,nod);

Snod=reshape(F.S(MUA.connectivity,1),Nele,nod);
Bnod=reshape(F.B(MUA.connectivity,1),Nele,nod);

unod=reshape(F.ub(MUA.connectivity,1),Nele,nod);
vnod=reshape(F.vb(MUA.connectivity,1),Nele,nod);

Psi_xnod=reshape(Psi_x(MUA.connectivity,1),Nele,nod);
Psi_ynod=reshape(Psi_y(MUA.connectivity,1),Nele,nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),Nele,nod);
nnod=reshape(F.n(MUA.connectivity,1),Nele,nod);

Cnod=reshape(F.C(MUA.connectivity,1),Nele,nod);
mnod=reshape(F.m(MUA.connectivity,1),Nele,nod);

rhonod=reshape(F.rho(MUA.connectivity,1),Nele,nod);

Fuu=zeros(Nele,nod,nod);
Fvv=zeros(Nele,nod,nod);
Fuv=zeros(Nele,nod,nod);

% Int loop and first node loop fully vectorized 
for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);

    Dx=reshape(Deriv(:,1,:),Nele,nod);
    Dy=reshape(Deriv(:,2,:),Nele,nod);

    h=hnod*fun;

    B=Bnod*fun;
    S=Snod*fun;
    H=S-B;

    u=unod*fun;
    v=vnod*fun;

    % values of the adjoint fields at the integration point
    Psi_x_int=Psi_xnod*fun;
    Psi_y_int=Psi_ynod*fun;

    n=nnod*fun;
    A=AGlennod*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;

    m=mnod*fun;
    C=Cnod*fun;
    C(C<CtrlVar.Cmin)=CtrlVar.Cmin;

    rho=rhonod*fun;

    hf=F.rhow*H./rho;
    He=HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
    G=He;

    % --- gradients at the integration point -----------------------------

    dudx=sum(Dx.*unod,2);
    dudy=sum(Dy.*unod,2);
    dvdx=sum(Dx.*vnod,2);
    dvdy=sum(Dy.*vnod,2);

    dPsi_xdx=sum(Dx.*Psi_xnod,2);
    dPsi_xdy=sum(Dy.*Psi_xnod,2);
    dPsi_ydx=sum(Dx.*Psi_ynod,2);
    dPsi_ydy=sum(Dy.*Psi_ynod,2);

    % --- constitutive quantities, all Nele x 1 --------------------------

    E_1=4.*dudx+2.*dvdy;
    E_2=dudy+dvdx;
    F_1=4.*dvdy+2.*dudx;

    e_eff=sqrt(dudx.^2+dvdy.^2+dudx.*dvdy+0.25.*(dvdx+dudy).^2+Eps0.^2);
    u_e=sqrt(u.^2+v.^2+u0.^2);

    % one vector-exponent power for the viscosity family
    tEta=(e_eff./A).^(1./n);              % = A^(-1/n) e_eff^(1/n)
    p3=tEta./e_eff.^3;                    % = A^(-1/n) e_eff^((1-3n)/n)
    p5=p3./e_eff.^2;                      % = A^(-1/n) e_eff^((1-5n)/n)

    c1=(1-n)./(8.*n).*p3;
    c2=((1-n).*(1-3.*n))./(32.*n.^2).*p5;

    % one vector-exponent power for the sliding family
    tBeta=(u_e./(C+C0)).^(1./m);          % = (C+C0)^(-1/m) u_e^(1/m)
    q3=tBeta./u_e.^3;                     % = (C+C0)^(-1/m) u_e^((1-3m)/m)
    q5=q3./u_e.^2;                        % = (C+C0)^(-1/m) u_e^((1-5m)/m)

    b1=(1-m)./m.*q3;
    b2=(1-m)./m.*q5;

    % --- contraction weights --------------------------------------------

    Weta =h.*(E_1.*dPsi_xdx+E_2.*(dPsi_xdy+dPsi_ydx)+F_1.*dPsi_ydy);
    Wbeta=G.*(u.*Psi_x_int+v.*Psi_y_int);

    Q1=4.*dPsi_xdx+2.*dPsi_ydy;
    Q2=dPsi_xdy+dPsi_ydx;
    Q3=2.*dPsi_xdx+4.*dPsi_ydy;

    detJw=detJ*MUA.weights(Iint);

    % quadrature weight folded in here, once, rather than onto Nele x nod arrays
    WC2=Weta.*c2.*detJw;
    WC1=Weta.*c1.*detJw;
    HC1=h.*c1.*detJw;

    % the complete sliding contribution collapses onto three Nele x 1 fields,
    % because the phi_i phi_k dependence is a pure scalar product
    Nuu=(Wbeta.*(b2.*((1-2.*m)./m.*u.^2+v.^2+u0.^2))+2.*G.*b1.*u.*Psi_x_int).*detJw;
    Nvv=(Wbeta.*(b2.*((1-2.*m)./m.*v.^2+u.^2+u0.^2))+2.*G.*b1.*v.*Psi_y_int).*detJw;
    Nuv=(Wbeta.*(((1-3.*m)./m).*b2.*u.*v)+G.*b1.*(v.*Psi_x_int+u.*Psi_y_int)).*detJw;

    % --- per-node arrays, Nele x nod ------------------------------------

    T =E_1.*Dx+E_2.*Dy;
    Sq=F_1.*Dy+E_2.*Dx;
    Ru=Q1.*Dx+Q2.*Dy;
    Rv=Q3.*Dy+Q2.*Dx;

    % --- node loop: Inod only, Knod vectorised --------------------------

    for Inod=1:nod

        Ti =T(:,Inod);
        Si =Sq(:,Inod);
        Rui=Ru(:,Inod);
        Rvi=Rv(:,Inod);
        dxi=Dx(:,Inod);
        dyi=Dy(:,Inod);

        PP=fun(Inod)*fun.';               % 1 x nod

        Fuu(:,:,Inod)=Fuu(:,:,Inod) ...
            + WC2.*(Ti.*T) + 4.*WC1.*(dxi.*Dx+0.25.*dyi.*Dy) ...
            + HC1.*(Ti.*Ru+T.*Rui) + Nuu.*PP;

        Fvv(:,:,Inod)=Fvv(:,:,Inod) ...
            + WC2.*(Si.*Sq) + 4.*WC1.*(dyi.*Dy+0.25.*dxi.*Dx) ...
            + HC1.*(Si.*Rv+Sq.*Rvi) + Nvv.*PP;

        Fuv(:,:,Inod)=Fuv(:,:,Inod) ...
            + WC2.*(Ti.*Sq) + WC1.*(2.*dxi.*Dy+dyi.*Dx) ...
            + HC1.*(Ti.*Rv+Sq.*Rui) + Nuv.*PP;

    end

end

%
%
% KFqq= [ Fuu Fuv ]
%       [ Fvu Fvv ]
%
% with F_vu(i,k) = F_uv(k,i), and the storage convention F**(:,Knod,Inod).
%
%

nnz4=4*nod*nod*Nele;
Iind=zeros(nnz4,1,'uint32');
Jind=zeros(nnz4,1,'uint32');
Xval=zeros(nnz4,1);

istak=0;
Nnodes=MUA.Nnodes;

for Inod=1:nod
    for Knod=1:nod

        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod);
        Xval(istak+1:istak+Nele)=Fuu(:,Knod,Inod);
        istak=istak+Nele;

        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod)+Nnodes;
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod)+Nnodes;
        Xval(istak+1:istak+Nele)=Fvv(:,Knod,Inod);
        istak=istak+Nele;

        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod)+Nnodes;
        Xval(istak+1:istak+Nele)=Fuv(:,Knod,Inod);
        istak=istak+Nele;

        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod)+Nnodes;
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod);
        Xval(istak+1:istak+Nele)=Fuv(:,Inod,Knod);   % F_vu(i,k) = F_uv(k,i)
        istak=istak+Nele;

    end
end

KFqq=sparse(Iind,Jind,Xval,2*Nnodes,2*Nnodes);

% symmetric by construction; a safeguard only
KFqq=(KFqq'+KFqq)/2;

if CtrlVar.Inverse.TestDirectAdjoint.isTrue

    FiniteDifferencesTestAndPlots(F,CtrlVar,MUA,BCs,Psi_x,Psi_y,KFqq,Nnodes);

end


end





function FiniteDifferencesTestAndPlots(F,CtrlVar,MUA,BCs,Psi_x,Psi_y,KFqq,Nnodes)
% note: These second-order derivatives are all identically equal to zero for n=1 and m=1

if all(F.m==1) && all(F.n==1)

    fprintf("Fqq: no real need to test these F_qq entries as they are all identically equal to zero for n=1 and m=1. \n")

end

CtrlVar.uvMatrixAssembly.Ronly = false;
CtrlVar.uvAssembly.ZeroFields  = false;
iColumn=randi(MUA.Nnodes);

hstep = 1e-6*max(sqrt(F.ub.*F.ub+F.vb.*F.vb)) ;

%% uu and vu
u0 = F.ub;

F.ub = u0; F.ub(iColumn) = F.ub(iColumn) - hstep;
[~,Kuv_minus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_minus = Kuv_minus.' * [Psi_x; Psi_y];

F.ub = u0; F.ub(iColumn) = F.ub(iColumn) + hstep;
[~,Kuv_plus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_plus = Kuv_plus.' * [Psi_x; Psi_y];

F.ub = u0;

Fqq_col_FD = (b_plus - b_minus)/(2*hstep);

Diff=norm(KFqq(:,iColumn) - Fqq_col_FD)/(norm(Fqq_col_FD)+eps);
fprintf("Fqq_v2: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

FiniteDifferenceApproximation=Fqq_col_FD;
KFqqColumn=full(KFqq(:,iColumn));

%% uv and vv
v0 = F.vb;

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) - hstep;
[~,Kuv_minus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_minus = Kuv_minus.' * [Psi_x; Psi_y];

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) + hstep;
[~,Kuv_plus] = uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
b_plus = Kuv_plus.' * [Psi_x; Psi_y];

F.vb = v0;

Fqq_col_FD = (b_plus - b_minus)/(2*hstep);

Diff=norm(KFqq(:,iColumn+Nnodes) - Fqq_col_FD)/(norm(Fqq_col_FD)+eps);
fprintf("Fqq_v2: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn+Nnodes,Diff)

FiniteDifferenceApproximation=[FiniteDifferenceApproximation;Fqq_col_FD];

KFqqColumn=[KFqqColumn;full(KFqq(:,iColumn+Nnodes))];

%% Plots
FqqTest=FindOrCreateFigure("Test: Fqq") ;

plot(KFqqColumn,FiniteDifferenceApproximation,"or") ; axis equal ;
hold on ;
plot([min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],[min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],"--k")

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")

title("$\mathcal{F}^{qq}_{uu,ik}$, $\mathcal{F}^{qq}_{uv,ik}$ , $\mathcal{F}^{qq}_{vu,ik}$ and $\mathcal{F}^{vv}_{vu,ik}$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

%%
end
