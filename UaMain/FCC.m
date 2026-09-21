function KFCC=FCC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y)

narginchk(7,7)

%%
%
%
% $$
% \mathcal{F}^{pp}_{\hat{C}\hat{C},lm} = \langle \Psi_x \mid \delta^2_{\hat{C}\hat{C}} \mathcal{F}_x[\phi_l,\phi_m] \rangle
%                          + \langle \Psi_y \mid \delta^2_{\hat{C}\hat{C}} \mathcal{F}_y[\phi_l,\phi_m] \rangle
% $$
%
% See FCC.m for the full documentation of the method and the log10
% conversion.
%
%
%%

ndim=2;
nNodes=MUA.Nnodes ;

Nele=MUA.Nele;
nod=MUA.nod;
nip=MUA.nip;

C0=CtrlVar.Czero;
u0=CtrlVar.SpeedZero;

h_node=reshape(F.h(MUA.connectivity,1),Nele,nod);   % Nele x nod
B_node=reshape(F.B(MUA.connectivity,1),Nele,nod);
S_node=reshape(F.S(MUA.connectivity,1),Nele,nod);
rho_node=reshape(F.rho(MUA.connectivity,1),Nele,nod);
hf_node=F.rhow*(S_node-B_node)./rho_node;

C_node=reshape(F.C(MUA.connectivity,1),Nele,nod);
m_node=reshape(F.m(MUA.connectivity,1),Nele,nod);

u_node=reshape(F.ub(MUA.connectivity,1),Nele,nod);
v_node=reshape(F.vb(MUA.connectivity,1),Nele,nod);

Psi_x_node=reshape(Psi_x(MUA.connectivity,1),Nele,nod);
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),Nele,nod);

W=zeros(Nele,nip);          % integration-point weights
P=zeros(nip,nod*nod);       % outer products of the shape functions

for Iint=1:nip

    fun=shape_fun(Iint,ndim,nod,MUA.points) ;
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

    % single vector-exponent power:
    %   (C+C0)^(-1/m-2) u_e^(1/m-1) = (u_e/(C+C0))^(1/m) / ( (C+C0)^2 u_e )
    CC=C+C0;
    tBeta=(u_e./CC).^(1./m);

    d2fdxdx=G.*((1+m)./m.^2).*tBeta./(CC.^2.*u_e).*(Psix.*u+Psiy.*v);

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

KFCC=sparse(Iind,Jind,Xval,nNodes,nNodes);

%%
% This function calculates:
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

    b_C=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);  % this function has already done the lin to log conversion
    D = spdiags(C_nodal*ln10, 0, nNodes, nNodes);
    KFCC = D*KFCC*D  + spdiags(b_C(:)*ln10, 0, nNodes, nNodes);

end


%% Test against finite differences

if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(MUA,CtrlVar,F,BCs,BCsAdjoint,Psi_x,Psi_y,KFCC);
end


end




function FiniteDifferenceTestAndPlots(MUA,CtrlVar,F,BCs,BCsAdjoint,Psi_x,Psi_y,KFCC)

iColumn=randi(MUA.Nnodes);  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.

CtrlVar.log10Derivatives=true;

C0=F.C;
logC0=log10(C0);

CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize=0.0001;

F.C=C0;

if CtrlVar.log10Derivatives
    dlogC=1e-5;
else
    deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(C0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;
    dC=deltaStep(iColumn);
end

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
fprintf("FCC_v2: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)

FCCTest=FindOrCreateFigure("Test: FCC") ; plot(KFCC(:,iColumn),HFCC_col_FD,"or") ; axis equal ;
hold on ;

Upper=max([HFCC_col_FD;KFCC(:,iColumn)]);
Lower=min([HFCC_col_FD;KFCC(:,iColumn)]);

plot([Lower Upper],[Lower Upper],"--k")

ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$\mathcal{F}^{pp}_{CC,lm} = \langle \Psi_x,\delta^2_{CC}\mathcal{F}_x[\phi_l,\phi_m] \rangle  + \langle \Psi_y , \delta^2_{CC}\mathcal{F}_y[\phi_l,\phi_m] \rangle$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")

drawnow

end
