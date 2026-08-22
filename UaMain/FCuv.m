




function [KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%%
%
% Calculates
%
%
% $$
% \mathcal{F}^{pq}_{Cu,li} =  \langle \Psi_x , \delta^2_{Cu}\mathcal{F}_x[\phi_l,\phi_i] \rangle   + \langle \Psi_y,\delta^2_{Cu}\mathcal{F}_y[\phi_l,\phi_i] \rangle
% $$
%
% $$ 
% \mathcal{F}^{pq}_{Cv,li} = \langle \Psi_x , \delta^2_{Cv}\mathcal{F}_x[\phi_l,\phi_i] \rangle  + \langle \Psi_y,\delta^2_{Cv}\mathcal{F}_y[\phi_l,\phi_i] \rangle
% $$
%
%
%
% For example using Weertman sliding law:
%
% $$ F^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi_i \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, u | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% $$ F^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi_i \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, v | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
% where
%
% $$\beta^2= (C+C_0)^{-1/m} \, (u^2+v^2+v_0^2)^{(1/m-1)/2} $$
%
% the relevant part is
%
% $$F_x= \mathcal{G} \beta^2 \, u $$
%
% $$F_y= \mathcal{G} \beta^2 \, v $$
%
%
% $$ \beta^2 = (C+C_0)^{-1/m} \, (u^2+v^2 + v_0)^{(1-m)/2m} $$
%
% $$F_x= \mathcal{G} \; (C+C_0)^{-1/m} \;  (u^2+v^2+v_0^2)^{(1-m)/2m}  \, u $$
%
% we have
%
% $$\delta_C F_x = -\frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; (u^2+v^2 + v_0^2)^{(1-m)/2m} \, u \; \delta C$$
%
% and
%
% $$\delta^2_{u \,C } F_x = -\frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; \left  (\frac{1-m}{m} (u^2+v^2 + v_0^2)^{(1-3m)/2m} \, u \, u +  (u^2+v^2 + v_0)^{(1-m)/2m} \right )    \; \delta C \, \delta u$$
%
%
%
% $$ \langle \Psi_x | \delta^2_{u C} F_x \rangle =
%  -\int  \Psi_x \,  \frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; \left  (\frac{1-m}{m} (u^2+v^2 + v_0)^{(1-3m)/2m} \, u \, u +  (u^2+v^2 + v_0^2)^{(1-m)/2m} \right )    \; \phi_i \, \phi_j \; dx \, dy
% $$
%
% which is a $n \times n$ matrix, which then needs to be multiplied with the $n \times n$ sensitivity matrix $\partial
% u/\partial C$
%
% There are two components to $F=(F^x,F^y)$ and two components to $q=(u,v)$, so we get
%
% $$ \langle \Psi_x | \delta^2_{u C} F_x \rangle \, \delta_c u + \langle  \Psi_x | \delta^2_{v C} F_x   \rangle \, \delta_C v
% +  \langle \Psi_y | \delta^2_{u C} F_y \rangle \, \delta_C u + \langle  \Psi_y | \delta^2_{v C} F_y   \rangle \, \delta_C v $$
%
%
%%


if ~contains(CtrlVar.Inverse.InvertFor,"logC",IgnoreCase=true)
    KFCu=[];
    KFCv=[];
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

C0=CtrlVar.Cmin;
u0=CtrlVar.SpeedZero;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod

Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hfnod=F.rhow*(Snod-Bnod)./rhonod;

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psixnod=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psiynod=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);



Hu=zeros(MUA.Nele,MUA.nod,MUA.nod);
Hv=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    % Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    h=hnod*fun;
    hf=hfnod*fun;

    u=unod*fun;
    v=vnod*fun;



    C=Cnod*fun;
    C(C<C0)=C0;

    m=mnod*fun;

    Psix=Psixnod*fun;
    Psiy=Psiynod*fun;

    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

    detJw=detJ*MUA.weights(Iint);

    % dcFx=G.*(-1./m).* (C+C0).^(-1./m-1) .*  (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)).*u ;
    % dcFy=G.*(-1./m).* (C+C0).^(-1./m-1) .*  (u.^2+v.^2 + v0^2).^((1-m)./(2.*m)).*v ;


    ducFx=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + u0^2).^((1-m)./(2.*m)-1) .*u.*u + (u.^2+v.^2 + u0^2).^((1-m)./(2.*m))  ) ;
    dvcFx=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + u0^2).^((1-m)./(2.*m)-1) .*v.*u  ) ;

    ducFy=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + u0^2).^((1-m)./(2.*m)-1) .*u.*v  ) ;
    dvcFy=G.*(-1./m).* (C+C0).^(-1./m-1) .* ((1./m-1) .* (u.^2+v.^2 + u0^2).^((1-m)./(2.*m)-1) .*v.*v + (u.^2+v.^2 + u0^2).^((1-m)./(2.*m))  ) ;


    l_d2FdudC=Psix.*ducFx+Psiy.*ducFy;  % raw kernel long log(C) scaling here at all
    l_d2FdvdC=Psix.*dvcFx+Psiy.*dvcFy;

    %% Testing
    % u_e=sqrt(u.*u + v.*v+u0^2) ;
    %
    %
    % Pi = -((1-m)*log(10))./(m.^2) .* C .* (C+C0).^(-1./m-1) .* u_e.^((1-3.*m)./m) ;
    % Sigma = -(log(10))./m .* C .* (C+C0).^(-1./m-1) .* u_e.^((1-m)./m) ;
    %
    % Ku=G.*((Pi.*u.*u+ Sigma).*Psix  +   Pi.*u.*v        .*Psiy) ;
    % Kv=G.*( Pi.*u.*v        .*Psix  +  (Pi.*v.*v+Sigma) .*Psiy ) ;
    % norm(Ku-l_d2FdudC) % is zero
    % norm(Kv-l_d2FdvdC) % is zero

    %%

    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            Hu(:,Inod,Jnod)=Hu(:,Inod,Jnod) + l_d2FdudC .*fun(Inod) .*fun(Jnod).*detJw;
            Hv(:,Inod,Jnod)=Hv(:,Inod,Jnod) + l_d2FdvdC .*fun(Inod) .*fun(Jnod).*detJw;

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

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        HuVal(istak+1:istak+MUA.Nele)=Hu(:,Inod,Jnod);
        HvVal(istak+1:istak+MUA.Nele)=Hv(:,Inod,Jnod);


        istak=istak+MUA.Nele;

    end
end

KFCu=sparse(Iind,Jind,HuVal,nNodes,nNodes) ;
KFCv=sparse(Iind,Jind,HvVal,nNodes,nNodes) ;

ln10 = log(10);
D_C = spdiags(F.C(:)*ln10, 0, nNodes, nNodes);

KFCu = D_C * KFCu;   % KFCu here is the RAW matrix from the assembly above
KFCv = D_C * KFCv;

KFCu = -KFCu;         % existing overall sign-convention fix, unchanged
KFCv = -KFCv;

%% Test Against Finite-Differences



if CtrlVar.Inverse.TestDirectAdjoint.isTrue

    FiniteDifferenceTestAndPlots(MUA,F,CtrlVar,BCs,BCsAdjoint,Psi_x,Psi_y,KFCu,KFCv);


end



end




function FiniteDifferenceTestAndPlots(MUA,F,CtrlVar,BCs,BCsAdjoint,Psi_x,Psi_y,KFCu,KFCv)


iColumn=randi(MUA.Nnodes); % just do the finite-difference comparison for this column of the Hessian contribution Fpp.



%% FCu test
u0 = F.ub;
v0=  F.vb;
hstep = 1e-6*max(abs(u0)+abs(v0));


F.ub = u0; F.ub(iColumn) = F.ub(iColumn) - hstep;
bMinus = dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

F.ub = u0; F.ub(iColumn) = F.ub(iColumn) + hstep;
bPlus = dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

F.ub = u0;

HCu_col_FD = (bPlus - bMinus)/(2*hstep);
Diff=norm(KFCu(:,iColumn) - HCu_col_FD)/norm(HCu_col_FD);
fprintf("FCu: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)


% FCv test
v0 = F.vb;


F.vb = v0; F.vb(iColumn) = F.vb(iColumn) - hstep;
bMinus = dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) + hstep;
bPlus = dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

F.vb = v0;

HCv_col_FD = (bPlus - bMinus)/(2*hstep);
Diff=norm(KFCv(:,iColumn) - HCv_col_FD)/norm(HCv_col_FD);
fprintf("FCv: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

%% Plots
FCuTest=FindOrCreateFigure("Test: FCu") ;

plot(KFCu(:,iColumn),HCu_col_FD,"or") ; axis equal ;
hold on ;
plot([min(HCu_col_FD) max(HCu_col_FD)],[min(HCu_col_FD) max(HCu_col_FD)],"--k")

% ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$\langle \Psi_x | \delta^2_{u C} F_x \rangle $",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")


FCvTest=FindOrCreateFigure("Test: FCv") ;

plot(KFCv(:,iColumn),HCv_col_FD,"or") ; axis equal ;
hold on ;
plot([min(HCv_col_FD) max(HCv_col_FD)],[min(HCv_col_FD) max(HCv_col_FD)],"--k")

% ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$\langle  \Psi_y | \delta^2_{v C} F_y   \rangle $",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

% 
% fprintf("FCuv: Inspect in debugger and then continue [F5]: \n")
% keyboard

%%
end