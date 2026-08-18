


function KFCC=FCC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%%
%
% Template for calculating terms such as
%
%
% $$
% \mathcal{F}^{pp}_{\hat{A}\hat{A},lm} = \langle \Psi_x \mid \delta^2_{\hat{A}\hat{A}} \mathcal{F}_x[\phi_l,\phi_m] \rangle
%                          + \langle \Psi_y \mid \delta^2_{\hat{A}\hat{A}} \mathcal{F}_y[\phi_l,\phi_m] \rangle
% $$
%
% and
%
%
% $$
% \mathcal{F}^{pp}_{\hat{C}\hat{C},lm} = \langle \Psi_x \mid \delta^2_{\hat{C}\hat{C}} \mathcal{F}_x[\phi_l,\phi_m] \rangle
%                          + \langle \Psi_y \mid \delta^2_{\hat{C}\hat{C}} \mathcal{F}_y[\phi_l,\phi_m] \rangle
% $$
%
% where
%
% $$ \hat{A}=\log_{10}(A) $$
%
% $$ \hat{C}=\log_{10}(C) $$
%
%
% The assembly is done with respect to $A$ and $C$ by expanding them as
%
% $$A(x,y) = \sum_i  A_i \, \phi(x,y) $$
% $$C(x,y) = \sum_i  C_i \, \phi(x,y) $$
%
%
% The FE assembly calculates these derivatives with respect to $A$ and $C$ and then afterwards the conversion to log space is
% done on the nodal vectors/matrices.
%
% Input variables:
%
% MUA  : a structure that contains data related the the FE mesh
%        For example:
%
% MUA.coordinates   : is as nNodes times 2 array with the x,y coordinates of the nodes
% MUA.connectivity  : contains the mesh connectivity, for a 3-nod element this would be a nEle times 3 array
% MUA.Nele          : Number of elements in the mesh
% MUA.Nnodes        : Number of nodes in the mesh
%
%
% F                 : a structure that contains all the field variables, for example:
% F.ub              : are the x velocity components (ub)
% F.vb              : the y velocity components (vb)
% F.h               : is the ice thickness (h)
%
%
%%



ndim=2;
nNodes=MUA.Nnodes ;


eta0=CtrlVar.etaZero;
Eps0=CtrlVar.EpsZero;
C0=CtrlVar.Czero;
u0=CtrlVar.SpeedZero;

h_node=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
B_node=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
S_node=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rho_node=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hf_node=F.rhow*(S_node-B_node)./rho_node;

C_node=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
m_node=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

A_node=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
n_node=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

u_node=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
v_node=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_x_node=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);



H=zeros(MUA.Nele,MUA.nod,MUA.nod);


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    % these are all values at integration points

    h=h_node*fun;
    hf=hf_node*fun;
    u=u_node*fun;     % u velocity at integration point Iint
    v=v_node*fun;     % v velocity at integration point
    u_e=sqrt(u.*u + v.*v+u0^2) ;

    C=C_node*fun;
    m=m_node*fun;



    Psix=Psi_x_node*fun;
    Psiy=Psi_y_node*fun;


    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % this is an external function


    dudx=zeros(MUA.Nele,1);
    dudy=zeros(MUA.Nele,1);

    dvdx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);

    dPsi_x_dx=zeros(MUA.Nele,1);
    dPsi_y_dx=zeros(MUA.Nele,1);
    dPsi_x_dy=zeros(MUA.Nele,1);
    dPsi_y_dy=zeros(MUA.Nele,1);

    % examples of how to calculate derivatives at integration points of various variables
    % not all of these may actually be needed in each case
    for Inode=1:MUA.nod

        dudx=dudx+Deriv(:,1,Inode).*u_node(:,Inode);
        dvdy=dvdy+Deriv(:,2,Inode).*v_node(:,Inode);
        dudy=dudy+Deriv(:,2,Inode).*u_node(:,Inode) ;
        dvdx=dvdx+Deriv(:,1,Inode).*v_node(:,Inode) ;

        dPsi_x_dx=dPsi_x_dx+Deriv(:,1,Inode).*Psi_x_node(:,Inode);
        dPsi_y_dx=dPsi_y_dx+Deriv(:,1,Inode).*Psi_y_node(:,Inode);
        dPsi_x_dy=dPsi_x_dy+Deriv(:,2,Inode).*Psi_x_node(:,Inode);
        dPsi_y_dy=dPsi_y_dy+Deriv(:,2,Inode).*Psi_y_node(:,Inode);

    end


    %   d2dfxdx = ...      % this is the point-wise expression at the integration point Iint

    d2fdxdx=G.* ((1+m)./m.^2) .* (C+C0).^(-1./m-2) .* u_e.^(1./m-1) .*(Psix.*u+Psiy.*v);


    detJw=detJ*MUA.weights(Iint);

    for Inode=1:MUA.nod  % This is a loop over nodes/basis vectors
        for Jnode=1:MUA.nod % This is another loop over nodes/basis vectors


            phi_i=fun(Inode);                 % \phi_i
            phi_j=fun(Jnode) ;                % \phi_j
            % dphidx_i=Deriv(:,1,Inode);        % d\phi_i/dx
            % dphidy_i=Deriv(:,2,Inode);        % d\phi_i/dy
            % dphidx_j=Deriv(:,1,Jnode);        % d\phi_j/dx
            % dphidy_j=Deriv(:,2,Jnode);        % d\phi_j/dy



            H(:,Inode,Jnode)=H(:,Inode,Jnode) + d2fdxdx .*phi_i .*phi_j.*detJw;

        end
    end
end


% Now we build the matrix, collect all values and indices into vectors and only make one call to the sparse function
Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);


istak=0;
for Inode=1:MUA.nod


    for Jnode=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inode);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnode);
        Xval(istak+1:istak+MUA.Nele)=H(:,Inode,Jnode);
        istak=istak+MUA.Nele;

    end
end


KFCC=sparse(Iind,Jind,Xval,nNodes,nNodes);

%%
% This function calculates:
%
%
% $$ b_l(C)=\langle  \delta_{C_i} F^x | \Psi_x \rangle + \langle  \delta_{C_i} F^y | \Psi_y \rangle $$
%
%%


% now I have the matrix Hpp, this includes second-order derivatives with respect to either A or C, but I now need those
% derivatives with respect log_{10}(A) and/or log_{10}(C)


if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    %% Convert Hpp (currently in C-space, since d2fdxdx above is the C-kernel)
    %  to log10(C)-space, using the nodal post-assembly formula:
    %
    %  H^{(Chat)} = D * Hpp * D  +  diag( b_C .* C .* ln10^2 )
    %
    %  where D = diag(C * ln10), applied entirely on the GLOBAL nodal vectors
    %  (F.C, b_C), not inside the integration-point loop.

    ln10 = log(10);
    C_nodal = F.C(:);

    b_C=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
    D = spdiags(C_nodal*ln10, 0, nNodes, nNodes);
    KFCC = D*KFCC*D  - spdiags(b_C(:)*ln10, 0, nNodes, nNodes);

end

KFCC=-KFCC;

%% Test against finite differences



if CtrlVar.Inverse.TestDirectAdjoint.isTrue

    iColumn=randi(MUA.Nnodes);  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.

    CtrlVar.log10Derivatives=true;

    %   KFCC=FCC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    C0=F.C;
    logC0=log10(C0);

    CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize=0.0001;

    % comparison

    F.C=C0;
    % perturbation

    if CtrlVar.log10Derivatives
        % perturb in log(C) space

        dlogC=1e-5;

    else
        % perturb in linear space
        deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(C0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;
        dC=deltaStep(iColumn);
    end

    % comparison

    if CtrlVar.log10Derivatives
        F.C(iColumn)=10^(logC0(iColumn)-dlogC);
    else
        F.C(iColumn)=F.C(iColumn)-dC;
    end

    FpMinus=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.C=C0;

    if CtrlVar.log10Derivatives
        F.C(iColumn)=10^(logC0(iColumn)+dlogC);
    else
        F.C(iColumn)=F.C(iColumn)+dC;
    end

    FpPlus=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    if CtrlVar.log10Derivatives
        HFCC_col_FD=(FpPlus-FpMinus)/(2*dlogC);  % Finite difference estimate of the iColumn of the Hessian

    else
        HFCC_col_FD=(FpPlus-FpMinus)/(2*dC);
    end



    Diff=norm(KFCC(:,iColumn)-HFCC_col_FD)/norm(HFCC_col_FD);
    fprintf("normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

    FCCTest=FindOrCreateFigure("FCC Test") ; plot(KFCC(:,iColumn),HFCC_col_FD,"or") ; axis equal ;
    hold on ;
    plot([min(HFCC_col_FD) max(HFCC_col_FD)],[min(HFCC_col_FD) max(HFCC_col_FD)],"--k")

    % ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

    xlabel("Direct-Adjoint",Interpreter="latex")  ;
    ylabel("Finite difference",Interpreter="latex")
    title("$\mathcal{F}^{pp}_{CC,lm} = \langle \Psi_x,\delta^2_{CC}\mathcal{F}_x[\phi_l,\phi_m] \rangle  + \langle \Psi_y , \delta^2_{CC}\mathcal{F}_y[\phi_l,\phi_m] \rangle$",Interpreter="latex")
    subtitle("Comparison is here for one random column",Interpreter="latex")

    drawnow
    input("Inspect, and then press RET to continue")


end



end





