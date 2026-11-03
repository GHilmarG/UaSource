function [f0,vs]=hElementKernel(MUA,Deriv,DetJ,conn,...
    h0v,h1v,u0v,v0v,u1v,v1v,rhov,a0v,a1v,dadhv,LSFMaskv,hBCsMaskv,...
    CtrlVar,Ronly,rPtr,rIdx,invPtr,invIdx)

%% Element-based h assembly kernel (mass conservation, velocities given), compiled with MATLAB Coder (31 Oct 2026); see README.md
%
% Called by hMatrixAssemblyCompiled.m through the MEX files built by BuildUaCompiled.m (entry points hKernelFull.m, hKernelRonly.m).
% Replaces the loop over integration points in MassContinuityEquationAssembly.m by a loop over elements (parfor, OpenMP), with the
% loop over integration points inside. The terms are those of MassContinuityEquationAssembly.m (the h row of the uvh system, with the
% velocities given), written for one element and one integration point; the pointwise functions (SUPGtau,
% LevelSetMethodMassBalanceFeedback, ThicknessPenaltyMassBalanceFeedback) are called with scalar arguments.
%
% Inputs:
%   MUA        : Nele, nod, nip, points, weights, EleAreas
%   Deriv, DetJ: MUA.Deriv and MUA.DetJ, passed as separate arrays (fields of a structure input are copied into a MEX file, arrays are not)
%   conn       : MUA.connectivity; the nodal values of each element are gathered from the nodal vectors (Nnodes x 1) here
%   h0v,h1v,u0v,v0v,u1v,v1v : thickness and velocities at the beginning (0) and end (1) of the time step
%   rhov       : density (F1.rho)
%   a0v,a1v    : mass balance (as+ab) at the beginning and end of the time step; dadhv its nodal derivative with respect to h
%                (mass-balance geometry feedback, otherwise zero)
%   LSFMaskv   : level-set mask (F1.LSFMask.NodesOut), used with CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback
%   hBCsMaskv  : nodes with thickness boundary conditions or active thickness constraints (no penalty in elements containing them)
%   CtrlVar    : run-time options and parameters (MakeKernelParamsFlex.m)
%   Ronly      : compile-time constant (residual only)
%   rPtr, rIdx : from CompiledAssemblyPattern.m, for summing the element residuals into the nodal residual
%   invPtr, invIdx : from CompiledAssemblyPattern.m, for summing the element-matrix values into the nonzeros
%
% Outputs:
%   f0 : the residual, Nnodes x 1 (K dh = -f0), summed in the kernel in the same order as accumarray(P.iR,...)
%   vs : the nonzero values of K, in CSC order (empty if Ronly)
%
% The lumped nodal thickness penalty (CtrlVar.ThicknessPenaltyNodal) and dFdt are not computed here (see MassContinuityEquationAssembly.m).
% Outputs that are completely overwritten are not zero-filled (coder.nullcopy). Interpreted, this function is slow.
%
%%

Nele=MUA.Nele; nod=MUA.nod; nip=MUA.nip; ndim=2;
assert(nod<=10) ; assert(nip<=37) ;   % upper bounds, so that the per-element arrays can be allocated statically

theta=CtrlVar.theta; dt=CtrlVar.dt;
beta0=CtrlVar.SUPG.beta0; SpeedZero=CtrlVar.SpeedZero;
tauOption=CtrlVar.uvh.SUPG.tauCode;   % numeric code, see SUPGtauCode.m and MakeKernelParamsFlex.m (the h assembly uses the uvh options)
tauMultiplier=CtrlVar.uvh.SUPG.tauMultiplier;

% options, constant over the loops
UseLSF=CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback;
UsePenalty=isfield(CtrlVar,'ThicknessPenalty') && CtrlVar.ThicknessPenalty && ~(isfield(CtrlVar,'ThicknessPenaltyNodal') && CtrlVar.ThicknessPenaltyNodal);

% shape functions at the integration points, the same for all elements
N=zeros(nod,nip);
for Iint=1:nip
    N(:,Iint)=shape_fun(Iint,ndim,nod,MUA.points);
end
weights=MUA.weights;

% element residuals and element-matrix values, element-major: Rmat(i,e), Xmat(L,e) with L = j + nod*(i-1)
Rmat=coder.nullcopy(zeros(nod,Nele));
if ~Ronly
    Xmat=coder.nullcopy(zeros(nod*nod,Nele));
else
    Xmat=zeros(nod*nod,0);
end

parfor e=1:Nele
    idx=conn(e,:);   % nodes of element e

    % nodal values of this element, 1 x nod
    h0n=h0v(idx).'; h1n=h1v(idx).'; u0n=u0v(idx).'; v0n=v0v(idx).'; u1n=u1v(idx).'; v1n=v1v(idx).'; rhon=rhov(idx).';
    a0n=a0v(idx).'; a1n=a1v(idx).'; dadhn=dadhv(idx).';
    ConstrainedElement=any(hBCsMaskv(idx).'>0);   % element contains a node with h BCs or an active thickness constraint
    l=sqrt(2*MUA.EleAreas(e));

    rh=zeros(1,nod);         % element residual
    kh=zeros(nod,nod);       % element matrix, kh(j,i): column node j, row node i

    for Iint=1:nip

        fun=N(:,Iint); funR=fun.';
        Dx=reshape(Deriv(e,1,:,Iint),1,nod);
        Dy=reshape(Deriv(e,2,:,Iint),1,nod);
        detJw=DetJ(e,Iint)*weights(Iint);

        % values at the integration point
        h0int=h0n*fun; h1int=h1n*fun;
        ub0int=u0n*fun; vb0int=v0n*fun;
        ub1int=u1n*fun; vb1int=v1n*fun;
        rhoint=rhon*fun;
        a0int=a0n*fun; a1int=a1n*fun;
        da1dhint=0;   % the derivative of the user-defined mass balance is a nodal quantity (dadhn)

        if UseLSF
            LM=LSFMaskv(idx).'*fun;
            [abLSF,dadhLSF]=LevelSetMethodMassBalanceFeedback(CtrlVar,LM,h1int);
            a1int=a1int+abLSF(1);
            da1dhint=da1dhint+dadhLSF(1);
        end

        if UsePenalty
            [aPenalty1,daPenaltydh1]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,h1int);
            a1int=a1int+aPenalty1(1)*(~ConstrainedElement);      % no penalty in elements containing h-constrained nodes
            da1dhint=da1dhint+daPenaltydh1(1)*(~ConstrainedElement);
        end

        % derivatives at the integration point
        dh1dx=Dx*h1n.'; dh1dy=Dy*h1n.';
        dh0dx=Dx*h0n.'; dh0dy=Dy*h0n.';
        exx0=Dx*u0n.'; eyy0=Dy*v0n.';
        drhodx=Dx*rhon.'; drhody=Dy*rhon.';
        exx1=Dx*u1n.'; eyy1=Dy*v1n.';

        speed0=sqrt(ub0int*ub0int+vb0int*vb0int+SpeedZero^2);
        tau=SUPGtau(CtrlVar,speed0,l,dt,tauOption,tauMultiplier);
        tauSUPGint=beta0*tau(1);

        q1xdx=rhoint*exx1*h1int+rhoint*ub1int*dh1dx+drhodx*ub1int*h1int;
        q1ydy=rhoint*eyy1*h1int+rhoint*vb1int*dh1dy+drhody*vb1int*h1int;
        q0xdx=rhoint*exx0*h0int+rhoint*ub0int*dh0dx+drhodx*ub0int*h0int;
        q0ydy=rhoint*eyy0*h0int+rhoint*vb0int*dh0dy+drhody*vb0int*h0int;

        SUPG=funR+0.5*tauSUPGint*(ub0int*Dx+vb0int*Dy);    % 1 x nod (i), as in the uvh assembly

        % residual: rho (h1-h0) + dt (theta div(q1) + (1-theta) div(q0)) - dt rho ((1-theta) a0 + theta a1)
        qterm=dt*(theta*q1xdx+(1-theta)*q0xdx+theta*q1ydy+(1-theta)*q0ydy);
        dhdt=rhoint*(h1int-h0int);
        accterm=dt*rhoint*((1-theta)*a0int+theta*a1int);
        rh=rh+((dhdt+qterm-accterm)*detJw)*SUPG;

        if ~Ronly
            % d(residual_i)/dh_j, 1 x nod (j)
            Kj=(rhoint*funR ...
                -dt*theta*rhoint*(da1dhint+dadhn).*funR ...
                +dt*theta*(rhoint*exx1*funR+drhodx*ub1int*funR+rhoint*ub1int*Dx ...
                          +rhoint*eyy1*funR+drhody*vb1int*funR+rhoint*vb1int*Dy))*detJw;
            for i=1:nod
                kh(:,i)=kh(:,i)+(SUPG(i)*Kj).';
            end
        end

    end

    Rmat(:,e)=rh.';
    if ~Ronly
        Xmat(:,e)=kh(:);   % element-major (contiguous for each element), L = j + nod*(i-1)
    end

end

% Sum the element residuals into the nodal residual, in parallel over the rows (each written by one thread), in the order given by rIdx
% (the same order as accumarray(P.iR,...), see CompiledAssemblyPattern.m).
neq=numel(rPtr)-1;
f0=coder.nullcopy(zeros(neq,1));
parfor r=1:neq
    s=0;
    for p=rPtr(r)+1:rPtr(r+1)
        s=s+Rmat(rIdx(p));
    end
    f0(r)=s;
end

% Sum the element-matrix values into the nonzeros of the sparse matrix, in parallel over the nonzeros, in the order given by invIdx.
if ~Ronly
    nk=numel(invPtr)-1;
    vs=coder.nullcopy(zeros(nk,1));
    parfor k=1:nk
        s=0;
        for p=invPtr(k)+1:invPtr(k+1)
            s=s+Xmat(invIdx(p));
        end
        vs(k)=s;
    end
else
    vs=zeros(0,1);
end
end
