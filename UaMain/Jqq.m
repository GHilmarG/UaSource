



function d2Jduvduv=Jqq(CtrlVar,MUA,F,BCs,Meas)

narginchk(5,5)

%%
%
% Calculates the $$J^{qq}$$ term, where $q$ is $u$ and $v$.
%
% There are number of terms in the cost function $J$, but the only ones which are explicit functions of $u$ and $v$ are:
%
% The $u$ misfit term:
%
% $$J_{u} = \frac{1}{2}\left(u_{\mathrm{obs}}-f_{u}\right)^{T}C_u^{-1}\left(u_{\mathrm{obs}}-f_v\right) $$
%
% The $v$ misfit term:
%
% $$J_{v} = \frac{1}{2}\left(v_{\mathrm{obs}}-f_{v}\right)^{T}C_v^{-1}\left(v_{\mathrm{obs}}-f_v\right) $$
%
% and the $\dot{h}$ misfit term:
%
% $$J_{\dot{h}}= \frac{1}{2\mathcal{A}} \int \left ( \dot{h}_{\mathrm{obs}} - a + \frac{1}{\rho} \left ( \partial_x (\rho h u)
% + \partial_y (\rho h v) \right ) \right )^2  \epsilon_{\dot{h}}^{-2}(x,y) \; dx \, dy
% $$
%
% Derivation and some useful quantities:
%
%
% $$
% r = \dot h_\mathrm{obs} - a + \left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)u + h\,\partial_xu + \left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)v + h\,\partial_yv
% $$
%
%
% $$
% \delta^2_{uu}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\Big(\delta_ur[\phi_i]\Big)\Big(\delta_ur[\phi_k]\Big)\epsilon_{\dot h}^{-2}\,dx\,dy
% + \underbrace{\frac{1}{\mathcal A}\int r\,\delta^2_{uu}r[\phi_i,\phi_k]\,\epsilon_{\dot h}^{-2}\,dx\,dy}_{=0,\ \mathrm{since}\;\delta^2_{uu}r\equiv0}
% $$
%
% $$
% \kappa_x := \partial_xh+\frac{h}{\rho}\partial_x\rho, \qquad \kappa_y := \partial_yh+\frac{h}{\rho}\partial_y\rho
% $$
%
% $$
% \delta_ur[\phi_i] = \kappa_x\phi_i+h\,\partial_x\phi_i, \qquad \delta_vr[\phi_i] = \kappa_y\phi_i+h\,\partial_y\phi_i
% $$
%
% $$
% \delta^2_{uu}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_x\phi_i+h\partial_x\phi_i\big)\big(\kappa_x\phi_k+h\partial_x\phi_k\big)\,dx\,dy
% $$
%
% $$
% \delta^2_{vv}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_y\phi_i+h\partial_y\phi_i\big)\big(\kappa_y\phi_k+h\partial_y\phi_k\big)\,dx\,dy
% $$
%
% $$
% \delta^2_{uv}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_x\phi_i+h\partial_x\phi_i\big)\big(\kappa_y\phi_k+h\partial_y\phi_k\big)\,dx\,dy
% $$
%
% $$
% J^{qq}_{\dot h} = \frac{1}{\mathcal A}\int \epsilon_{\dot h}^{-2}\;\big(\delta_qr\big)\otimes\big(\delta_qr\big)\;dx\,dy
% $$
%
% $$
% r = \dot h_\mathrm{obs} - a + \left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)u + h\,\partial_xu + \left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)v + h\,\partial_yv
% $$
%%



ndim=2;
Area=MUA.Area;


d2Jduvduv=0;


if contains(CtrlVar.Inverse.Measurements,'-uv-','IgnoreCase',true)

    %% calculate the contributions from the u and v misfit terms
    %
    % $$ 0.5 (u_{obs}-u)^T M (u_{obs}-u)/(err^2 Area) $$
    %
    % and
    %
    % $$ 0.5 (v_{obs}-u)^T M (v_{obs}-v)/(err^2 Area)$$
    %
    % Here M is the mass matrix, which is a part of the MUA structure, M=MUA.M
    %
    % the err in u and v are the diagonals of the covariance matrices Meas.usCov and Meas.vsCov
    %
    % This misfit term is always included, so I don't need to check if -uv- is in the string CtrlVar.Inverse.Measurements,
    % (although at some later stage I might consider allowing for only dh/dt measurements being used).
    %
    %
    uErr2=spdiags(Meas.usCov);
    d2Iduu=MUA.M./uErr2/Area;  % MUA.M is the mass matrix, OK this is not quite correct, only correct if errors are not spatially variable...

    vErr2=spdiags(Meas.vsCov);
    d2Idvv=MUA.M./vErr2/Area;
    d2Jduvduv=blkdiag(d2Iduu,d2Idvv);  % This is a block matrix, no mixing uv terms here



end


if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)

    %% OK, now we calculate the contribution due to dh/dt misfit term $$J_{\dot{h}}$$ shown above.
    %
    % This is only needed if we are using measurements of $\dot{h}$ in the inversion.
    %
    %

    hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
    %anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
    unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

    rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

    Juu=zeros(MUA.Nele,MUA.nod,MUA.nod);
    Jvu=zeros(MUA.Nele,MUA.nod,MUA.nod);
    Juv=zeros(MUA.Nele,MUA.nod,MUA.nod);
    Jvv=zeros(MUA.Nele,MUA.nod,MUA.nod);

    %dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dhdtErr=sqrt(spdiags(Meas.dhdtCov));
    dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),MUA.Nele,MUA.nod);


    for Iint=1:MUA.nip

        fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
        detJ=MUA.DetJ(:,Iint);
        Deriv=MUA.Deriv(:,:,:,Iint);

        hint=hnod*fun;

        % aint=anod*fun;
        %
        % uint=unod*fun;
        % vint=vnod*fun;



        dhdtErrint=dhdtErrnod*fun;
        % dhdtMeasint=dhdtMeasnod*fun;

        rhoint=rhonod*fun;



        dhdx=zeros(MUA.Nele,1);
        dhdy=zeros(MUA.Nele,1);


        drhodx=zeros(MUA.Nele,1);
        drhody=zeros(MUA.Nele,1);

        dudx=zeros(MUA.Nele,1);
        dudy=zeros(MUA.Nele,1);

        dvdx=zeros(MUA.Nele,1);
        dvdy=zeros(MUA.Nele,1);



        for Inod=1:MUA.nod

            dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
            dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);


            drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
            drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);

            dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
            dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
            dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod) ;
            dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod) ;


        end

        %   hdot=aint-(rhoint.*dhdx.*uint+rhoint.*hint.*dudx+drhodx.*hint.*uint+rhoint.*dhdy.*vint+rhoint.*hint.*dvdy+drhody.*hint.*vint)./rhoint ;
        %  r=(hdot-dhdtMeasint)./dhdtErrint;


        kx=dhdx+hint.*drhodx./rhoint;
        ky=dhdy+hint.*drhody./rhoint;

        detJw=detJ*MUA.weights(Iint);


        Multy=1./(Area.*dhdtErrint.^2);
        % denominator

        for Inod=1:MUA.nod
            for Knod=1:MUA.nod

                phi_k=fun(Knod);
                phi_i=fun(Inod) ;
                dphidx_i=Deriv(:,1,Inod);
                dphidy_i=Deriv(:,2,Inod);
                dphidx_k=Deriv(:,1,Knod);
                dphidy_k=Deriv(:,2,Knod);



                J_uu_ik = (kx.*phi_i+hint.*dphidx_i).*(kx.*phi_k+hint.*dphidx_k).*Multy;
                J_vv_ik = (ky.*phi_i+hint.*dphidy_i).*(ky.*phi_k+hint.*dphidy_k).*Multy;
                J_uv_ik = (kx.*phi_i+hint.*dphidx_i).*(ky.*phi_k+hint.*dphidy_k).*Multy;
                J_vu_ik = (ky.*phi_i+hint.*dphidy_i).*(kx.*phi_k+hint.*dphidx_k).*Multy;

                Juu(:,Inod,Knod)=Juu(:,Inod,Knod) + J_uu_ik.*detJw;
                Juv(:,Inod,Knod)=Juv(:,Inod,Knod) + J_uv_ik.*detJw;
                Jvu(:,Inod,Knod)=Jvu(:,Inod,Knod) + J_vu_ik.*detJw;
                Jvv(:,Inod,Knod)=Jvv(:,Inod,Knod) + J_vv_ik.*detJw;


            end
        end
    end


    Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
    Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
    Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

    istak=0;

    %
    %
    % KJqq= [ Juu Juv ]
    %       [ Jvu Jvv ]
    %
    %
    %
    %

    Nnodes=MUA.Nnodes;
    for Inod=1:MUA.nod
        for Knod=1:MUA.nod

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod);
            Xval(istak+1:istak+MUA.Nele)=Juu(:,Inod,Knod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+Nnodes;
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod)+Nnodes;
            Xval(istak+1:istak+MUA.Nele)=Jvv(:,Inod,Knod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod)+Nnodes;
            Xval(istak+1:istak+MUA.Nele)=Juv(:,Inod,Knod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+Nnodes;
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Knod);
            Xval(istak+1:istak+MUA.Nele)=Jvu(:,Inod,Knod);
            istak=istak+MUA.Nele;


        end
    end



    d2Jduvduv=sparse(Iind,Jind,Xval,2*Nnodes,2*Nnodes)+d2Jduvduv ; %% adding all up.

end


if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,Meas,d2Jduvduv)
end




end







function FiniteDifferenceTestAndPlots(CtrlVar,MUA,F,BCs,Meas,d2Jduvduv)


iColumn=randi(MUA.Nnodes);
hstep = 1e-7*max(sqrt(F.ub.*F.ub+F.vb.*F.vb));


%% uu uv
u0 = F.ub;


F.ub = u0; F.ub(iColumn) = F.ub(iColumn) - hstep;
[dJduvMinus]=Jq(CtrlVar,MUA,F,BCs,Meas);

F.ub = u0; F.ub(iColumn) = F.ub(iColumn) + hstep;
[dJduvPlus]=Jq(CtrlVar,MUA,F,BCs,Meas);


F.ub = u0;
d2Jduvduv_col_FD = (dJduvPlus - dJduvMinus)/(2*hstep);   % length 2*Nnodes: top half -> column of F^{qq}_{uu}, bottom half -> column of F^{qq}_{vu}


FiniteDifferenceApproximationA=d2Jduvduv_col_FD;
d2JduvduvColumnA=full(d2Jduvduv(:,iColumn));

Diff=norm(d2JduvduvColumnA - FiniteDifferenceApproximationA)/(norm(d2JduvduvColumnA)+eps);
fprintf("Jqq: normalized norm of difference for uu/uv column %i is %g \n",iColumn,Diff)

%% vu vv
v0 = F.vb;



F.vb = v0; F.vb(iColumn) = F.vb(iColumn) - hstep;
[dJduvMinus]=Jq(CtrlVar,MUA,F,BCs,Meas);

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) + hstep;
[dJduvPlus]=Jq(CtrlVar,MUA,F,BCs,Meas);


F.vb = v0;
d2Jduvduv_col_FD = (dJduvPlus - dJduvMinus)/(2*hstep);   % length 2*Nnodes: top half -> column of F^{qq}_{uu}, bottom half -> column of F^{qq}_{vu}


FiniteDifferenceApproximationB=d2Jduvduv_col_FD;
d2JduvduvColumnB=full(d2Jduvduv(:,iColumn+MUA.Nnodes));

Diff=norm(d2JduvduvColumnB - FiniteDifferenceApproximationB)/(norm(d2JduvduvColumnB)+eps);
fprintf("Jqq: normalized norm of difference for uu/uv column %i is %g \n",iColumn,Diff)


%% comparison 
FiniteDifferenceApproximation=[FiniteDifferenceApproximationA;FiniteDifferenceApproximationB];
d2JduvduvColumn=[d2JduvduvColumnA;d2JduvduvColumnB];
%% Plots
JqqTest=FindOrCreateFigure("Test: Jqq") ;

plot(d2JduvduvColumn,FiniteDifferenceApproximation,"or") ; axis equal ;
hold on ;
plot([min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],[min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],"--k")

% ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("$\delta^2_{uu,uv,vu,vv} J$",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")


title("$J^{qq}_{uu,ik}$, $J^{qq}_{uv,ik}$ , $J^{qq}_{vu,ik}$ and $J^{vv}_{vu,ik}$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

% fprintf("Fqq: Inspect in debugger and then continue: [F5] \n")
% keyboard

%%
end





