


function KFAA=FAA(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

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

% 
% eta0=CtrlVar.etaZero;
% Eps0=CtrlVar.EpsZero;
% C0=CtrlVar.Czero;
% u0=CtrlVar.SpeedZero;

h_node=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod


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

    n=n_node*fun;
    A=A_node*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;


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

    exx=dudx;
    eyy=dvdy;
    exy=0.5*(dudy+dvdx);

     [~,~,~,~,d2etadAdA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy);


     %eEff=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));
     %d2etadAdA=real(  A.^(-1./n-2) .* eEff.^(1./n-1) .* 0.5.*(1./(n.^2)+(1./n)) );

     Temp=h.*( (4*exx+2*eyy).*dPsi_x_dx + 2*exy.*dPsi_x_dy + (4*eyy+2*exx).*dPsi_y_dy + 2*exy.*dPsi_y_dx); 
     
     d2fdxdx=d2etadAdA.*Temp;  % this is the point-wise expression at the integration point Iint, may have to be edited further


    detJw=detJ*MUA.weights(Iint);

    for Inode=1:MUA.nod  % This is a loop over nodes/basis vectors
        for Jnode=1:MUA.nod % This is another loop over nodes/basis vectors


            phi_i=fun(Inode);                 % \phi_i
            phi_j=fun(Jnode) ;                % \phi_j
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


KFAA=sparse(Iind,Jind,Xval,nNodes,nNodes);


% now I have the matrix Hpp, this includes second-order derivatives with respect to either A or C, but I now need those
% derivatives with respect log_{10}(A) and/or log_{10}(C)

if CtrlVar.log10Derivatives

    %% Convert Hpp (currently in A-space, since d2fdxdx above is the A-kernel)
   
    ln10 = log(10);
    A_nodal = F.AGlen(:);
    b_A=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
    D = spdiags(A_nodal*ln10, 0, nNodes, nNodes);
    KFAA = -D*KFAA*D  + spdiags(b_A(:)*ln10, 0, nNodes, nNodes);


end

%KFAA=-KFAA;

%% Test against finite differences


TestAgainstFiniteDifferences=true;  % This is a flag to test against finite-differences, get-rid of this in production

if TestAgainstFiniteDifferences

  
    iColumn=5;  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.

    CtrlVar.log10Derivatives=true;

 %   KFAA=FAA(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    A0=F.AGlen;
    logA0=log10(A0);



    CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize=0.0001;

    % comparison

    F.AGlen=A0;
    % perturbation

    if CtrlVar.log10Derivatives
        % perturb in log(A) space

        dlogA=1e-6;

    else
        % perturb in linear space
        deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(A0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;
        dA=deltaStep(iColumn);
    end

    % comparison

    if CtrlVar.log10Derivatives
        F.AGlen(iColumn)=10^(logA0(iColumn)-dlogA);
    else
        F.AGlen(iColumn)=F.AGlen(iColumn)-dA;
    end

    FpMinus=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    F.AGlen=A0;

    if CtrlVar.log10Derivatives
        F.AGlen(iColumn)=10^(logA0(iColumn)+dlogA);
    else
        F.AGlen(iColumn)=F.AGlen(iColumn)+dA;
    end

    FpPlus=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

    if CtrlVar.log10Derivatives
        HFAA_col_FD=(FpPlus-FpMinus)/(2*dlogA);  % Finite difference estimate of the iColumn of the Hessian

    else
        HFAA_col_FD=(FpPlus-FpMinus)/(2*dA);
    end



    Diff=norm(KFAA(:,iColumn)-HFAA_col_FD)/norm(HFAA_col_FD);
    fprintf("normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

    FAATest=FindOrCreateFigure("FAA Test") ; plot(KFAA(:,iColumn),HFAA_col_FD,"or") ; axis equal ;
    hold on ;
    plot([min(HFAA_col_FD) max(HFAA_col_FD)],[min(HFAA_col_FD) max(HFAA_col_FD)],"--k")

    % ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

    xlabel("Direct-Adjoint",Interpreter="latex")  ;
    ylabel("Finite difference",Interpreter="latex")
    title("$\mathcal{F}^{pp}_{AA,lm} = \langle \Psi_x,\delta^2_{AA}\mathcal{F}_x[\phi_l,\phi_m] \rangle  + \langle \Psi_y , \delta^2_{AA}\mathcal{F}_y[\phi_l,\phi_m] \rangle$",Interpreter="latex")
    subtitle("Comparison is here for one random column",Interpreter="latex")

    fprintf("                                                  AD                        FD \n")
    l=10 ; [(1:l)' full(KFAA(1:l,iColumn)) HFAA_col_FD(1:l) full(KFAA(1:l,iColumn))-HFAA_col_FD(1:l)]
  

  
end













end





