function KFAA=FAA(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%%
%
%
% $$
% \mathcal{F}^{pp}_{\hat{A}\hat{A},lm} = \langle \Psi_x \mid \delta^2_{\hat{A}\hat{A}} \mathcal{F}_x[\phi_l,\phi_m] \rangle
%                          + \langle \Psi_y \mid \delta^2_{\hat{A}\hat{A}} \mathcal{F}_y[\phi_l,\phi_m] \rangle
% $$
%
% See FAA.m for the full documentation of the method and the log10
% conversion.
%
% The element contribution is
%
%    H(:,i,j) = sum_Iint  w_Iint  fun_Iint(i) fun_Iint(j)
%
% with w = d2fdxdx.*detJw an Nele x 1 field and fun(i) a scalar.  The
% element matrix is therefore a mass matrix weighted by w, i.e. a rank-one
% outer product at each integration point.  
%%

ndim=2;
nNodes=MUA.Nnodes ;

Nele=MUA.Nele;
nod=MUA.nod;
nip=MUA.nip;

h_node=reshape(F.h(MUA.connectivity,1),Nele,nod);   % Nele x nod

A_node=reshape(F.AGlen(MUA.connectivity,1),Nele,nod);
n_node=reshape(F.n(MUA.connectivity,1),Nele,nod);

u_node=reshape(F.ub(MUA.connectivity,1),Nele,nod);
v_node=reshape(F.vb(MUA.connectivity,1),Nele,nod);

Psi_x_node=reshape(Psi_x(MUA.connectivity,1),Nele,nod);
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),Nele,nod);

W=zeros(Nele,nip);          % integration-point weights
P=zeros(nip,nod*nod);       % outer products of the shape functions

for Iint=1:nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    Dx=reshape(Deriv(:,1,:),Nele,nod);
    Dy=reshape(Deriv(:,2,:),Nele,nod);

    % these are all values at integration points

    h=h_node*fun;

    n=n_node*fun;
    A=A_node*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;

    dudx=sum(Dx.*u_node,2);
    dudy=sum(Dy.*u_node,2);
    dvdx=sum(Dx.*v_node,2);
    dvdy=sum(Dy.*v_node,2);

    dPsi_x_dx=sum(Dx.*Psi_x_node,2);
    dPsi_x_dy=sum(Dy.*Psi_x_node,2);
    dPsi_y_dx=sum(Dx.*Psi_y_node,2);
    dPsi_y_dy=sum(Dy.*Psi_y_node,2);

    exx=dudx;
    eyy=dvdy;
    exy=0.5*(dudy+dvdx);

    [~,~,~,~,d2etadAdA]=EffectiveViscositySSTREAM(CtrlVar,A,n,exx,eyy,exy);

    Temp=h.*( (4*exx+2*eyy).*dPsi_x_dx + 2*exy.*dPsi_x_dy + (4*eyy+2*exx).*dPsi_y_dy + 2*exy.*dPsi_y_dx);

    d2fdxdx=d2etadAdA.*Temp;  % point-wise expression at the integration point Iint

    detJw=detJ*MUA.weights(Iint);

    W(:,Iint)=d2fdxdx.*detJw;
    P(Iint,:)=reshape(fun*fun.',1,nod*nod);

end

H=reshape(W*P,Nele,nod,nod);


% Now we build the matrix, collect all values and indices into vectors and only make one call to the sparse function
Iind=zeros(nod*nod*Nele,1,'uint32');
Jind=zeros(nod*nod*Nele,1,'uint32');
Xval=zeros(nod*nod*Nele,1);

istak=0;
for Inode=1:nod
    for Jnode=1:nod

        Iind(istak+1:istak+Nele)=MUA.connectivity(:,Inode);
        Jind(istak+1:istak+Nele)=MUA.connectivity(:,Jnode);
        Xval(istak+1:istak+Nele)=H(:,Inode,Jnode);
        istak=istak+Nele;

    end
end

KFAA=sparse(Iind,Jind,Xval,nNodes,nNodes);


% now I have the matrix Hpp, this includes second-order derivatives with respect to either A or C, but I now need those
% derivatives with respect log_{10}(A) and/or log_{10}(C)

if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    %% Convert Hpp (currently in A-space, since d2fdxdx above is the A-kernel)

    ln10 = log(10);
    A_nodal = F.AGlen(:);
    b_A=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
    D = spdiags(A_nodal*ln10, 0, nNodes, nNodes);
    KFAA = D*KFAA*D  + spdiags(b_A(:)*ln10, 0, nNodes, nNodes);

end


%% Test against finite differences

if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(MUA,CtrlVar,F,BCs,BCsAdjoint,Psi_x,Psi_y,KFAA);
end


end




function FiniteDifferenceTestAndPlots(MUA,CtrlVar,F,BCs,BCsAdjoint,Psi_x,Psi_y,KFAA)

iColumn=randi(MUA.Nnodes);  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.

CtrlVar.log10Derivatives=true;

A0=F.AGlen;
logA0=log10(A0);

CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize=0.0001;

F.AGlen=A0;

if CtrlVar.log10Derivatives
    dlogA=1e-6;
else
    deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(A0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;
    dA=deltaStep(iColumn);
end

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
fprintf("FAA_v2: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

FAATest=FindOrCreateFigure("Test: FAA") ; plot(KFAA(:,iColumn),HFAA_col_FD,"or") ; axis equal ;
hold on ;
plot([min(HFAA_col_FD) max(HFAA_col_FD)],[min(HFAA_col_FD) max(HFAA_col_FD)],"--k")

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$\mathcal{F}^{pp}_{AA,lm} = \langle \Psi_x,\delta^2_{AA}\mathcal{F}_x[\phi_l,\phi_m] \rangle  + \langle \Psi_y , \delta^2_{AA}\mathcal{F}_y[\phi_l,\phi_m] \rangle$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

end
