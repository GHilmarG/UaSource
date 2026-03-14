

function [Kdh1du,Kdh1dv]=dhduv(CtrlVar,MUA,F,BCs)

narginchk(4,4)



%%
%
% $$
% F_{h_1}(u(A)) = h_1(u(A))- h_0 +\Delta t \, ( \partial_x (u(A) h_0) - a ) = 0
% $$
%
% and
%
% $$
% 0=\frac{\partial F_{h_1}}{\partial h_1} \frac{\partial h_1}{\partial u} + \frac{\partial F_{1}}{\partial u}
% $$
%
%
% or
%
% $$
%  \langle h_1 \, | \, \phi_k  \rangle  - \langle h_0 \, | \, \phi_k \rangle  + \Delta t \, \langle  \partial_x (u h_0) - a
%  \, | \, \phi_k \rangle  =0
% $$
%
% and
%
% $$
%  \langle \phi_i \, | \, \phi_k  \rangle  \; \delta_{u} h_1   + \Delta t \, \langle  \partial_x (\phi_i h_0) \, | \, \phi_k \rangle  =0
% $$
%
%
%
%% Calculates the sensitivity matrix dh/duv
%
%
% The $k$-column contains the response in u and v to a perturbation in $A_k$
%
%
% $$\left[\begin{array}{cccc}
% \partial h^1_1 /\partial u_1  & \partial h^1_1 /\partial u_2  & \ldots & \partial h^1_1 /\partial u_n  \\
% \partial h^1_2 /\partial u_1  & \partial h^1_2 /\partial u_2  & \ldots & \partial h^1_2 /\partial u_n  \\
%              .              &              .              &  .  &    .                          \\
% \partial h^1_1 /\partial u_n  & \partial h^1_2 /\partial h^1_n & \ldots & \partial h^1_n /\partial u_n  \\
% \end{array}\right] $$
%
%
% Boundary conditions:
%
% I assume that where we have boundary conditions on u, v and h, the corresponding sensitives are zero, as any change in p
% (A, B, C) can not affect q (u,v, dh/dt) at those nodes.
%
% This boundary conditions are linear, and the problem itself is linear. So this can be solved as
%
% $$
% \left ( \begin{array}{cc}
%       \frac{\partial F}{\partial q}  & L^T \\
%       L       &  0
% \end{array} \right )
% \left ( \begin{array}{c}
%       \xi  \\
%        0
% \end{array} \right )
% =
% - \left ( \begin{array}{c}
%       \frac{\partial F}{\partial p}  \\
%        0
% \end{array} \right )
% $$
%
%
% see also: duvdAFunc.m, duvdBFunc.m, duvdCFunc.m, dFuvdA.m, dFuvdB.m, dFuvdC.m, TestSensitivityMatrixCalculations.m
%
%%



ndim=2;
nNodes=MUA.Nnodes ;
dt=F.dt;


hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

dFh1du=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFh1dv=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFh1h1=zeros(MUA.Nele,MUA.nod,MUA.nod);

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    h=hnod*fun;

    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dudx=dudx+Deriv(:,1,Inod).*ubnod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vbnod(:,Inod);

        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);


    end

    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            dFh1du(:,Inod,Jnod)=dFh1du(:,Inod,Jnod) + dt*(h.* Deriv(:,1,Jnod) + fun(Jnod).*dhdx ).*fun(Inod).*detJw;
            dFh1dv(:,Inod,Jnod)=dFh1dv(:,Inod,Jnod) + dt*(h.* Deriv(:,2,Jnod) + fun(Jnod).*dhdy ).*fun(Inod).*detJw;
            dFh1h1(:,Inod,Jnod)=dFh1h1(:,Inod,Jnod) + fun(Jnod).*fun(Inod) .*detJw;  % I don't really need to do this as this is the Mass Matrix

        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');

dFhdotduXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
dFhdotdvXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
dFhdotdhdotXval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);

        dFhdotduXval(istak+1:istak+MUA.Nele)=dFh1du(:,Inod,Jnod);
        dFhdotdvXval(istak+1:istak+MUA.Nele)=dFh1dv(:,Inod,Jnod);
        dFhdotdhdotXval(istak+1:istak+MUA.Nele)=dFh1h1(:,Inod,Jnod);

        istak=istak+MUA.Nele;

    end
end


dFh1du=sparseUA(Iind,Jind,dFhdotduXval,nNodes,nNodes);
dFh1dv=sparseUA(Iind,Jind,dFhdotdvXval,nNodes,nNodes);
dFh1h1=sparseUA(Iind,Jind,dFhdotdhdotXval,nNodes,nNodes);




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


hRhs=repmat(hRhs,1,size(dFh1du,2));
x0=[] ; y0=hRhs*0;
%dFh1h1=MUA.M; 
% I could do this in one solve by expanding the right-hand side, missing the minus sign
Kdh1du=solveKApeSymmetric(-dFh1h1,hL,dFh1du,hRhs,x0,y0,CtrlVar);
Kdh1dv=solveKApeSymmetric(-dFh1h1,hL,dFh1dv,hRhs,x0,y0,CtrlVar);

%%

TestSensitivityMatrixCalculations=true;

if TestSensitivityMatrixCalculations

    NodeTest=1045;
    NodeTest=1200;
    UserVar=[];



    dh1du=Kdh1du(:,NodeTest);
    dh1dv=Kdh1dv(:,NodeTest);

    u0=F.ub;
    v0=F.vb;
    du=u0(NodeTest)*0.01;
    dv=v0(NodeTest)*0.01;


    up=u0;
    up(NodeTest)=up(NodeTest)+du;
    F.ub=up; F.vb=v0;
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
    h1p=dhdt*F.dt ;

    um=u0;
    um(NodeTest)=um(NodeTest)-du;
    F.ub=um; F.vb=v0 ;
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
    h1m=dhdt*F.dt;

    dh1duPert=(h1p-h1m)/(2*du) ;

    
    vp=v0;
    vp(NodeTest)=vp(NodeTest)+dv;
    F.ub=u0;  F.vb=vp;
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
    h1p=dhdt*F.dt ;

    vm=v0;
    vm(NodeTest)=vm(NodeTest)-dv;
    F.ub=u0 ; F.vb=vm;
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
    h1m=dhdt*F.dt;

    dh1dvPert=(h1p-h1m)/(2*dv) ;


    % dh/du
    fighu=FindOrCreateFigure("dh/du comparision"); clf(fighu)
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dh1du,CreateNewFigure=false)  ; title("$dh/dv$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dh1duPert,CreateNewFigure=false) ; title("$dh/du$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dh1du-dh1duPert,CreateNewFigure=false) ; title("$dh/du$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";

    fighv=FindOrCreateFigure("dh/dv comparision"); clf(fighv)
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dh1dv,CreateNewFigure=false)  ; title("$dh/dv$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dh1dvPert,CreateNewFigure=false) ; title("$dh/dv$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dh1dv-dh1dvPert,CreateNewFigure=false) ; title("$dh/dv$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";


    fighhugrad=FindOrCreateFigure("dh/du gradient test") ;  clf(fighhugrad)
    plot(dh1du,dh1duPert,"or") ;
    hold on
    axis equal
    AX=axis;
    plot([min(dh1du) max(dh1du)],[min(dh1du) max(dh1du)],"--k") ;
    axis equal tight ;
    xlabel(" $dh/du$",Interpreter="latex")  ;
    ylabel("Finite difference $dh/du$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight ; box off
    title("Comparision between adjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')

    fighhugrad=FindOrCreateFigure("dh/dv gradient test") ;  clf(fighhugrad)
    plot(dh1dv,dh1dvPert,"or") ;
    hold on
    axis equal
    plot([min(dh1dv) max(dh1dv)],[min(dh1dv) max(dh1dv)],"--k") ;
    axis equal tight ;
    xlabel(" $dh/dv$",Interpreter="latex")  ;
    ylabel("Finite difference $dh/dv$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight ; box off
    title("Comparision between adjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')





end


end


















