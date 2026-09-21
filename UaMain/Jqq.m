function d2Jduvduv=Jqq(CtrlVar,MUA,F,BCs,Meas)

narginchk(5,5)

%%
%
% Refactored version of Jqq.m.  Mathematically identical.
%
% The $\dot h$ misfit contribution has the rank-one structure
%
%    delta_u r[phi_i] = kx phi_i + h dphidx_i  =:  Gx_i
%    delta_v r[phi_i] = ky phi_i + h dphidy_i  =:  Gy_i
%
%    J_uu(i,k) = Gx_i Gx_k W        W := detJw / (Area eps_hdot^2)
%    J_vv(i,k) = Gy_i Gy_k W
%    J_uv(i,k) = Gx_i Gy_k W
%    J_vu(i,k) = Gy_i Gx_k W       = J_uv(k,i)
%
% so per integration point only the two Nele x nod arrays Gx and Gy are
% needed.  The original rebuilt (kx phi_i + h dphidx_i) and its partner
% four times inside every (Inod,Knod) pair, i.e. 8*nod^2 times per
% integration point where 2*nod suffice.
%
% Changes relative to Jqq.m:
%
%  1) Gx, Gy hoisted out of the double node loop; the loop now runs over
%     Inod only, with Knod vectorised.
%
%  2) J_vu is not formed.  J_vu(i,k) = J_uv(k,i) exactly, and J_uu, J_vv
%     are symmetric in (i,k), so d2Jduvduv is symmetric by construction,
%     as it must be for a second variation.
%
%  3) The quadrature weight is folded into Multy once per integration
%     point rather than multiplied onto four Nele x nod arrays.
%
%  4) dudx, dudy, dvdx and dvdy were accumulated in the nodal loop but
%     never used - the only references are in commented-out lines.  They
%     are removed.  The nodal derivative accumulations that are used are
%     replaced by sum(...,2).
%
%  5) The index arrays of the assembly stage are preallocated to their
%     true length 4*nod*nod*Nele.  The original preallocated
%     nod*nod*Nele, i.e. a factor four short, so those arrays were being
%     grown in place.
%
% Storage convention: the element arrays are held as J**(:,Knod,Inod),
% i.e. [Nele , k , i], which makes the writes in the Inod loop
% contiguous.  The assembly stage indexes them accordingly.
%
% For the documentation of the underlying mathematics see Jqq.m.
%
%%

ndim=2;
Area=MUA.Area;

Nele=MUA.Nele;
nod=MUA.nod;

d2Jduvduv=0;


if contains(CtrlVar.Inverse.Measurements,'-uv-','IgnoreCase',true)

    %% contributions from the u and v misfit terms
    %
    % $$ 0.5 (u_{obs}-u)^T M (u_{obs}-u)/(err^2 Area) $$
    %
    % and
    %
    % $$ 0.5 (v_{obs}-u)^T M (v_{obs}-v)/(err^2 Area)$$
    %
    uErr2=spdiags(Meas.usCov);
    d2Iduu=MUA.M./uErr2/Area;  % MUA.M is the mass matrix, OK this is not quite correct, only correct if errors are not spatially variable...

    vErr2=spdiags(Meas.vsCov);
    d2Idvv=MUA.M./vErr2/Area;
    d2Jduvduv=blkdiag(d2Iduu,d2Idvv);  % This is a block matrix, no mixing uv terms here

end


if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)

    %% contribution due to the dh/dt misfit term $$J_{\dot{h}}$$

    hnod=reshape(F.h(MUA.connectivity,1),Nele,nod);
    rhonod=reshape(F.rho(MUA.connectivity,1),Nele,nod);

    Juu=zeros(Nele,nod,nod);
    Juv=zeros(Nele,nod,nod);
    Jvv=zeros(Nele,nod,nod);

    dhdtErr=sqrt(spdiags(Meas.dhdtCov));
    dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),Nele,nod);

    for Iint=1:MUA.nip

        fun=shape_fun(Iint,ndim,nod,MUA.points) ;
        detJ=MUA.DetJ(:,Iint);
        Deriv=MUA.Deriv(:,:,:,Iint);

        Dx=reshape(Deriv(:,1,:),Nele,nod);
        Dy=reshape(Deriv(:,2,:),Nele,nod);

        hint=hnod*fun;
        rhoint=rhonod*fun;
        dhdtErrint=dhdtErrnod*fun;

        dhdx=sum(Dx.*hnod,2);
        dhdy=sum(Dy.*hnod,2);

        drhodx=sum(Dx.*rhonod,2);
        drhody=sum(Dy.*rhonod,2);

        kx=dhdx+hint.*drhodx./rhoint;
        ky=dhdy+hint.*drhody./rhoint;

        detJw=detJ*MUA.weights(Iint);

        % quadrature weight folded in here, once
        W=detJw./(Area.*dhdtErrint.^2);

        % variation of the residual with respect to the nodal velocities,
        % Nele x nod
        Gx=kx.*fun.'+hint.*Dx;
        Gy=ky.*fun.'+hint.*Dy;

        for Inod=1:nod

            WGx=W.*Gx(:,Inod);
            WGy=W.*Gy(:,Inod);

            Juu(:,:,Inod)=Juu(:,:,Inod) + WGx.*Gx;
            Jvv(:,:,Inod)=Jvv(:,:,Inod) + WGy.*Gy;
            Juv(:,:,Inod)=Juv(:,:,Inod) + WGx.*Gy;

        end

    end

    %
    %
    % KJqq= [ Juu Juv ]
    %       [ Jvu Jvv ]
    %
    % with J_vu(i,k) = J_uv(k,i), and the storage convention J**(:,Knod,Inod).
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
            Xval(istak+1:istak+Nele)=Juu(:,Knod,Inod);
            istak=istak+Nele;

            Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod)+Nnodes;
            Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod)+Nnodes;
            Xval(istak+1:istak+Nele)=Jvv(:,Knod,Inod);
            istak=istak+Nele;

            Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod);
            Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod)+Nnodes;
            Xval(istak+1:istak+Nele)=Juv(:,Knod,Inod);
            istak=istak+Nele;

            Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inod)+Nnodes;
            Jind(istak+1:istak+Nele)=MUA.connectivity(:,Knod);
            Xval(istak+1:istak+Nele)=Juv(:,Inod,Knod);   % J_vu(i,k) = J_uv(k,i)
            istak=istak+Nele;

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
d2Jduvduv_col_FD = (dJduvPlus - dJduvMinus)/(2*hstep);

FiniteDifferenceApproximationA=d2Jduvduv_col_FD;
d2JduvduvColumnA=full(d2Jduvduv(:,iColumn));

Diff=norm(d2JduvduvColumnA - FiniteDifferenceApproximationA)/(norm(d2JduvduvColumnA)+eps);
fprintf("Jqq_v2: normalized norm of difference for uu/uv column %i is %g \n",iColumn,Diff)

%% vu vv
v0 = F.vb;

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) - hstep;
[dJduvMinus]=Jq(CtrlVar,MUA,F,BCs,Meas);

F.vb = v0; F.vb(iColumn) = F.vb(iColumn) + hstep;
[dJduvPlus]=Jq(CtrlVar,MUA,F,BCs,Meas);

F.vb = v0;
d2Jduvduv_col_FD = (dJduvPlus - dJduvMinus)/(2*hstep);

FiniteDifferenceApproximationB=d2Jduvduv_col_FD;
d2JduvduvColumnB=full(d2Jduvduv(:,iColumn+MUA.Nnodes));

Diff=norm(d2JduvduvColumnB - FiniteDifferenceApproximationB)/(norm(d2JduvduvColumnB)+eps);
fprintf("Jqq_v2: normalized norm of difference for vu/vv column %i is %g \n",iColumn+MUA.Nnodes,Diff)


%% comparison
FiniteDifferenceApproximation=[FiniteDifferenceApproximationA;FiniteDifferenceApproximationB];
d2JduvduvColumn=[d2JduvduvColumnA;d2JduvduvColumnB];

%% Plots
JqqTest=FindOrCreateFigure("Test: Jqq") ;

plot(d2JduvduvColumn,FiniteDifferenceApproximation,"or") ; axis equal ;
hold on ;
plot([min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],[min(FiniteDifferenceApproximation) max(FiniteDifferenceApproximation)],"--k")

xlabel("$\delta^2_{uu,uv,vu,vv} J$",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")

title("$J^{qq}_{uu,ik}$, $J^{qq}_{uv,ik}$ , $J^{qq}_{vu,ik}$ and $J^{vv}_{vu,ik}$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

%%
end
