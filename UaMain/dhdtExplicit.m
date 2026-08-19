
function [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs)
%%
% Calculates dh/dt from flux divergence as
%
% $$ \rho \, \dot{h}  = \rho \, a -  \nabla \cdot \mathbf{q} $$
%
% or
%
% $$ \dot{h}  = a -  \frac{1}{\rho} \nabla \cdot \mathbf{q} $$
%
% where
%
% $$\nabla \cdot \mathbf{q} = \partial_x (\rho h u ) + \partial_y (\rho h v) $$
%
%   [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F)
%
% uses u=F.ub, and hence only correct for plug flow, e.g. SSA
%
% Projects the values directly onto nodes.
%
% see also : ProjectFintOntoNodes
%
% Note: Homogenized  thickness (h) boundary conditions are applied to the dh/dt solve. So if thickness is set to some
% prescribed value at some nodes, dh/dt at those nodes will be forced to be equal to zero. And if thickness links/ties are defined
% between nodes, those same times will be applied to dh/dt.
%
% $$
% F_{\dot{h}} =\rho \, \dot{h}+\partial_x (\rho u h) + \partial_y (\rho v h) - \rho \, a = 0
% $$
%
% Solve for $\dot{h}$ as:
%
% $$ \delta_{\dot{h}} F_h + F_h(h_0) = 0 $$
%
% leading to
%
% $$ \partial_{\dot{h}} F_{\dot{h}} \, \Delta h = -F(h_0) $$
%
% $$
% \langle \rho \phi_i \, | \phi_j \rangle \, \Delta \dot{h}_i   =  \langle \rho \, a \, | \phi_j \rangle
% - \langle \partial_x (\rho u h ) \, | \phi_j \rangle
% - \langle \partial_y (\rho v h ) \, | \phi_j \rangle
% $$
%
%
% The NR-system to solve
%
% $$
% \langle \phi_i \, | \phi_k \rangle \, \Delta \dot{h}_i = -   \langle F_{\dot{h}} \, | \, \phi_j \rangle
% $$
%
% where
%
% $$ \dot{h} \to \dot{h} + \Delta \dot{h} $$
%
% until right-hand side small enough. However, as the equation is linear, in effect $\dot{h}$ it is solved directly, or in
% just one iteration.
%
%%


narginchk(5,5)

if nargin<5
    BCs=[];
end


ndim=2; dof=1; neq=dof*MUA.Nnodes;

anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


% [points,weights]=sample('triangle',MUA.nip,ndim);

b=zeros(MUA.Nele,MUA.nod);


if isempty(MUA.Deriv)  || anynan(MUA.Deriv)
    CtrlVar.CalcMUA_Derivatives=true;
    MUA=UpdateMUA(CtrlVar,MUA);
end

% vector over all elements for each integration point
for Iint=1:MUA.nip


    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);

    aint=anod*fun;
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    rhoint=rhonod*fun;

    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    drhodx=zeros(MUA.Nele,1);
    drhody=zeros(MUA.Nele,1);
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod

        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);

        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);

        drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
        drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);

    end

    % % these are the flux gradients at integration this integration point, Iint, for all elements
    % dqxdx=rhoint.*dudx.*hint+rhoint.*uint.*dhdx+drhodx.*uint.*hint;
    % dqydy=rhoint.*dvdy.*hint+rhoint.*vint.*dhdy+drhody.*vint.*hint;
    %
    % rhs=aint.*rhoint-dqxdx-dqydy;

    % these are the flux gradients at integration this integration point, Iint, for all elements, divided by rho
    dqxdx=dudx.*hint+uint.*dhdx+drhodx.*uint.*hint./rhoint;
    dqydy=dvdy.*hint+vint.*dhdy+drhody.*vint.*hint./rhoint;

    rhs=aint-dqxdx-dqydy;  % here I have already divided by rho



    detJw=detJ*MUA.weights(Iint);

    for Inod=1:MUA.nod

        b(:,Inod)=b(:,Inod)+rhs.*fun(Inod).*detJw;

    end



    % assemble right-hand side

    rh=sparseUA(neq,1);
    for Inod=1:MUA.nod
        rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b(:,Inod),neq,1);
    end

    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end

    %% BCs
    %
    % Since this is an explicit estimate, it seems right to use the BCs for h.
    % If h is fixed, the the explicit estimate should be dh/dt=0 for those nodes
    %
    % Also, if there is a nodal link (tie) for h, then the same nodal tie should be used for dh/dt
    %
    % So for all BCs.hFixed nodes, set BCs.hFixedValue=0, and use the h ties.


    [hL,hRhs]=createLh(MUA.Nnodes,BCs.hFixedNode,BCs.hFixedValue*0,BCs.hTiedNodeA,BCs.hTiedNodeB);

    %% Solve system
    %CtrlVar.SymmSolver='AugmentedLagrangian';


    if isempty(MUA.M)
        MUA.M=MassMatrix2D1dof(MUA);
    end

    
    x0=[] ; y0=hRhs*0;
    [dhdt,dhdtlambda]=solveKApeSymmetric(MUA.M,hL,rh,hRhs,x0,y0,CtrlVar);
    dhdt=full(dhdt);

    % Now there is an issue here regarding what to do about dhdt<0 when h<=thickmin

    I=(F.h<=CtrlVar.ThickMin)  & (dhdt< 0) ;
    dhdt(I)=0 ;


end