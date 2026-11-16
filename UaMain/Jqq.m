function d2Jduvduv=Jqq(CtrlVar,MUA,F,BCs,Meas)

narginchk(5,5)

persistent uvKey uvHessian   % (16 Nov 2026) cache of the (state-independent) velocity-misfit Hessian, see below

%%
%
% Refactored version of Jqq.m.  Mathematically identical.
%
% The $\dot h$ misfit contribution has the rank-one structure
%
%    delta_u r[phi_i] = kx phi_i + h dphidx_i  =:  Gx_i
%    delta_v r[phi_i] = ky phi_i + h dphidy_i  =:  Gy_i
%
%    J_uu(i,k) = Gx_i Gx_k W        W := detJw w / Area,  w = 1/eps_hdot^2 interpolated from the nodes
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
    % $$ I_{uv} = \frac{1}{2 \, Area} \left ( r_u^T M_{w_u} \, r_u + r_v^T M_{w_v} \, r_v \right ), \quad r_u=u-u_{obs}, \quad r_v=v-v_{obs} $$
    %
    % where M_w is the weighted mass matrix, (M_w)_{pq} = < w \phi_p , \phi_q >, with the nodal weights w=1/err^2 interpolated to
    % the integration points (see Misfit.m and Jq.m). This is the consistent discretisation of the misfit integral also when
    % the errors vary in space. The gradient is M_{w_u} r_u / Area, and the Hessian is therefore
    %
    % $$ M_{w_u} / Area $$
    %
    % and equivalently for v.
    %
    % (16 Nov 2026) This Hessian does not depend on the state (the misfit is quadratic in u and v): it only depends on the mesh, the
    % observation operator and the errors. It is therefore calculated once, and kept (persistent) for subsequent calls with the same
    % mesh and measurements, identified by a key calculated from these (UvHessianKey below). Nothing is stored in restart files: in a
    % restart run, or after "clear Jqq", it is calculated again at the first call. In the -points- formulation, with millions of
    % data, this was most of the time needed to assemble the Hessian blocks.
    Key=UvHessianKey(CtrlVar,MUA,Meas) ;

    if ~isempty(uvKey) && isequal(Key,uvKey)

        d2Iduu=[] ; d2Idvv=[] ;   % not needed, uvHessian is used below

    elseif is_uv_points_misfit(CtrlVar)

        % "-points-" formulation, the Hessian is O' W O, see PointDataMisfit.m (14 November 2026)
        [~,~,d2Iduu]=PointDataMisfit(Meas.uvO,Meas.usData,Meas.usErr,F.ub+F.ud);
        [~,~,d2Idvv]=PointDataMisfit(Meas.uvO,Meas.vsData,Meas.vsErr,F.vb+F.vd);

    else

        wu=1./full(diag(Meas.usCov));
        wv=1./full(diag(Meas.vsCov));
        [~,Muw]=GeneralizedMassAndStiffnessAssembly([],CtrlVar,MUA,wu,0,0,0);
        [~,Mvw]=GeneralizedMassAndStiffnessAssembly([],CtrlVar,MUA,wv,0,0,0);
        d2Iduu=(Muw+Muw.')/2/Area;     % symmetrised, as sparse() may sum duplicates in different order for (i,j) and (j,i)
        d2Idvv=(Mvw+Mvw.')/2/Area;

    end

    if isempty(uvKey) || ~isequal(Key,uvKey)
        uvHessian=blkdiag(d2Iduu,d2Idvv);  % This is a block matrix, no mixing uv terms here
        uvKey=Key ;
    end
    d2Jduvduv=uvHessian;

end


if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)

    %% contribution due to the dh/dt misfit term $$J_{\dot{h}}$$

    if dhdtMisfitMethod(CtrlVar)~="-nodal-integral-"

        % "-points-" and "-blocks-" formulations: Gauss-Newton term (A G)' W (A G) with G=[Gu Gv], see dhdtObservationMisfit.m.
        % There is no term proportional to the residuals as hdot is linear in u and v. (14 November 2026)
        [~,~,~,~,Mats]=dhdtObservationMisfit(CtrlVar,MUA,F,Meas,"Jacobians");
        AG=[Mats.AGu Mats.AGv];
        m=numel(Mats.w);
        d2Jduvduv=AG'*spdiags(Mats.w,0,m,m)*AG+d2Jduvduv;

    else

        hnod=reshape(F.h(MUA.connectivity,1),Nele,nod);
        rhonod=reshape(F.rho(MUA.connectivity,1),Nele,nod);

        Juu=zeros(Nele,nod,nod);
        Juv=zeros(Nele,nod,nod);
        Jvv=zeros(Nele,nod,nod);

        % data weights w=1/dhdtErr^2 at the nodes, interpolated to the integration points (as in EvaluateJhdotAndDerivatives.m)
        dhdtW=1./full(diag(Meas.dhdtCov));
        dhdtWnod=reshape(dhdtW(MUA.connectivity,1),Nele,nod);

        for Iint=1:MUA.nip

            fun=shape_fun(Iint,ndim,nod,MUA.points) ;
            detJ=MUA.DetJ(:,Iint);
            Deriv=MUA.Deriv(:,:,:,Iint);

            Dx=reshape(Deriv(:,1,:),Nele,nod);
            Dy=reshape(Deriv(:,2,:),Nele,nod);

            hint=hnod*fun;
            rhoint=rhonod*fun;
            dhdtWint=dhdtWnod*fun;

            dhdx=sum(Dx.*hnod,2);
            dhdy=sum(Dy.*hnod,2);

            drhodx=sum(Dx.*rhonod,2);
            drhody=sum(Dy.*rhonod,2);

            kx=dhdx+hint.*drhodx./rhoint;
            ky=dhdy+hint.*drhody./rhoint;

            detJw=detJ*MUA.weights(Iint);

            % quadrature weight folded in here, once
            W=detJw.*dhdtWint./Area;

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


function Key=UvHessianKey(CtrlVar,MUA,Meas)

% (16 Nov 2026) Key identifying the velocity-misfit Hessian: the formulation, the mesh, the observation operator (-points-) and the
% errors. Weighted sums are used so that the key changes when any of these change (including a permutation or a change of a single
% value), at the cost of a few vector operations and one sparse matrix-vector product.

N=MUA.Nnodes ;
xn=MUA.coordinates(:,1) ; yn=MUA.coordinates(:,2) ;
MeshPart=[N MUA.Nele MUA.nod sum(xn) sum(yn) sum(xn.*(1:N)') sum(yn.*(1:N)') sum(double(MUA.connectivity(:)).*(1:numel(MUA.connectivity))')] ;

if is_uv_points_misfit(CtrlVar)
    O=Meas.uvO ; m=size(O,1) ;
    Oi=O*(1:N)' ;                     % depends on the location of every non-zero of O
    wu=1./Meas.usErr(:).^2 ; wv=1./Meas.vsErr(:).^2 ;
    DataPart=[1 m nnz(O) sum(Oi) sum(Oi.*(1:m)') sum(wu) sum(wu.*(1:m)') sum(wv) sum(wv.*(1:m)')] ;
else
    wu=1./full(diag(Meas.usCov)) ; wv=1./full(diag(Meas.vsCov)) ;
    DataPart=[2 MUA.Area sum(wu) sum(wu.*(1:N)') sum(wv) sum(wv.*(1:N)')] ;
end

Key=[MeshPart DataPart] ;

end
