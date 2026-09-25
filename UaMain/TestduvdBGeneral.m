





function TestduvdBGeneral(CtrlVar,MUA,F,BCs,l,Nodes,dB,doPlots)

%% Tests the B sensitivity calculations against brute-force finite differences
%
% Tests
%
%   dFuvdBGeneral.m   and   duvdBGeneral.m
%
% for the general case where the ice may be grounded in some parts of the domain and afloat in others.
%
%   TestduvdBGeneral(CtrlVar,MUA,F,BCs,l)
%   TestduvdBGeneral(CtrlVar,MUA,F,BCs,l,Nodes,dB,doPlots)
%
%% What is, and is not, being tested
%
% duvdBGeneral.m returns the sensitivity (Jacobian) matrices
%
% $$ \left ( \frac{\partial u_i}{\partial B_j} \right ) , \qquad \left ( \frac{\partial v_i}{\partial B_j} \right ) $$
%
% Both indices are nodal indices. This is NOT the gradient of a scalar functional, and therefore no Riesz mapping, and
% no choice of metric, enters anywhere. A perturbation of the nodal value $B_j$ by $\delta$ gives, to first order,
% $\delta u_i = (\partial u_i / \partial B_j) \, \delta$, which is exactly what a finite difference of the forward
% solution measures. The comparison below is therefore like-for-like, column by column, with nothing to convert.
%
% (This is in contrast to dIdBqGeneral.m, which returns the components of the differential of the scalar cost
% function. Those are the components of a covector, and turning them into a search direction does require an inner
% product on the parameter space. That is what the Riesz mapping does, and why it must be switched off when testing
% that gradient against finite differences.)
%
% Note also that the columns of the Jacobian are not mesh-independent quantities: a unit perturbation of B at a single
% node is a different physical perturbation depending on the local element size. This does not affect the test, but it
% does mean that column magnitudes should not be compared across a non-uniform mesh without an area weighting.
%
%% The two tests performed
%
% # The contraction identity. dIdBqGeneral.m returns
%
% $$ \left ( \partial F^{uv} / \partial B \right )^{T} \left [ \Psi_x ; \Psi_y \right ] $$
%
% so for arbitrary $\Psi$ the matrix from dFuvdBGeneral.m, contracted with $\Psi$, must reproduce it to round-off.
% This requires no forward solves and tests every term of the sensitivity matrix at once.
%
% # Brute-force finite differences. For selected nodes, B is perturbed, the geometry is recalculated through
% Calc_bh_From_sBS.m, the non-linear uv system is re-solved, and a second-order centred difference is compared with
% the corresponding column of the Jacobian.
%
% By default one node is selected from each of the three grounding states (grounded, transition, afloat), chosen as
% the node with the largest velocity response within that state. Selecting nodes by index instead tends to pick nodes
% with a negligible response, for which the relative difference is a ratio of two numbers that are both effectively
% zero, and therefore meaningless.
%
%% Inputs
%
%   Nodes     optional. Nodes at which the finite-difference test is done. If empty, nodes are selected automatically
%             as described above.
%   dB        optional, default 0.01 m. The finite-difference step. Note that the grounding-line transition has a
%             width of order 1/CtrlVar.kH, and the step must be small compared with that: a step of 1 m gives
%             truncation errors of order 10 per cent at nodes near the grounding line.
%   doPlots   optional, default false.
%
%  see also: duvdBGeneral.m, dFuvdBGeneral.m, dIdBqGeneral.m, dGeometrydB.m
%
%%

narginchk(5,8)

if nargin<6 ; Nodes=[]   ; end
if nargin<7 || isempty(dB) ; dB=0.01 ; end
if nargin<8 || isempty(doPlots) ; doPlots=false ; end

CtrlVar.MapOldToNew.Test=false;
CtrlVar.Inverse.TestDirectAdjoint.isTrue=false;

%% make sure the geometry is the converged solution of the closure, then solve the forward problem

[F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow);

[~,~,F,l,KdFuvduv]= uv([],[],CtrlVar,MUA,BCs,F,l);

G=F.GF.node;

fprintf("\n TestduvdBGeneral: %i nodes.  grounded %4.1f%% , transition %4.1f%% , afloat %4.1f%% \n", ...
    MUA.Nnodes,100*mean(G>0.99),100*mean(G>=0.01 & G<=0.99),100*mean(G<0.01))

%% Test 1:  the contraction identity,  (dF/dB)' Psi  ==  dIdBqGeneral

K=dFuvdBGeneral(CtrlVar,MUA,F);

Psi_x=randn(MUA.Nnodes,1) ; Psi_y=randn(MUA.Nnodes,1) ;   % arbitrary adjoint fields

BCsAdjoint=BCs;   % not used by dIdBqGeneral other than being passed through
g=dIdBqGeneral(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

gK=K.'*[Psi_x;Psi_y];

fprintf("\n Test 1: (dFuvdB)'' Psi  against  dIdBqGeneral \n")
fprintf("   normalized difference = %g    (expect round-off) \n",norm(gK-g)/norm(g))

%% Test 2: finite differences of the forward solution

if isempty(Nodes)
    Nodes=iSelectNodes(CtrlVar,MUA,F,l,BCs,KdFuvduv,G);
end

[dudB,dvdB]=duvdBGeneral(CtrlVar,MUA,F,l,BCs,KdFuvduv,Nodes);

B0=F.B;

fprintf("\n Test 2: finite differences of the forward uv solution, step dB=%g m \n",dB)
fprintf("   node      GF        |dudB|       |FD u|      rel diff u      rel diff v \n")

for k=1:numel(Nodes)

    j=Nodes(k);

    Fp=F ; Fp.B(j)=B0(j)+dB ;
    [Fp.b,Fp.h,Fp.GF]=Calc_bh_From_sBS(CtrlVar,MUA,Fp.s,Fp.B,Fp.S,Fp.rho,Fp.rhow);
    [~,~,Fp,~]= uv([],[],CtrlVar,MUA,BCs,Fp,l);

    Fm=F ; Fm.B(j)=B0(j)-dB ;
    [Fm.b,Fm.h,Fm.GF]=Calc_bh_From_sBS(CtrlVar,MUA,Fm.s,Fm.B,Fm.S,Fm.rho,Fm.rhow);
    [~,~,Fm,~]= uv([],[],CtrlVar,MUA,BCs,Fm,l);

    dudBfd=(Fp.ub-Fm.ub)/(2*dB) ;
    dvdBfd=(Fp.vb-Fm.vb)/(2*dB) ;

    ru=norm(dudB(:,k)-dudBfd)/max(norm(dudBfd),realmin) ;
    rv=norm(dvdB(:,k)-dvdBfd)/max(norm(dvdBfd),realmin) ;

    fprintf("%7i  %7.4f  %11.4g  %11.4g  %14.4g  %14.4g",j,G(j),norm(dudB(:,k)),norm(dudBfd),ru,rv)

    if norm(dudBfd) < 1e-8*max(vecnorm(dudB))
        fprintf("   <- response negligible, ratio not meaningful")
    end
    fprintf("\n")

    if doPlots
        FigName=sprintf("dudB node %i",j);
        Fig=FindOrCreateFigure(FigName) ; clf(Fig)
        tiledlayout("flow")
        nexttile ; UaPlots(CtrlVar,MUA,F,[dudB(:,k) dvdB(:,k)],CreateNewFigure=false) ;   title("$(du/dB,dv/dB)$  model",Interpreter="latex") ; subtitle("")
        nexttile ; UaPlots(CtrlVar,MUA,F,[dudBfd,dvdBfd],CreateNewFigure=false) ;       title("$(du/dB,dv/dB)$  finite differences",Interpreter="latex") ; subtitle("")
        nexttile ; UaPlots(CtrlVar,MUA,F,[dudB(:,k)-dudBfd,dvdB(:,k)-dvdBfd],CreateNewFigure=false) ; title("difference") ;subtitle("")
    end

end

fprintf("\n")

end


function Nodes=iSelectNodes(CtrlVar,MUA,F,l,BCs,KdFuvduv,G)

% Select, from each grounding state, the node with the largest velocity response. A random subset of candidates is
% used so that only a modest number of columns of the Jacobian needs to be formed.

nSample=30 ;
Nodes=[] ;

Sets={ find(G>0.99) , find(G>=0.01 & G<=0.99) , find(G<0.01) } ;

for k=1:numel(Sets)

    I=Sets{k} ;
    if isempty(I) ; continue ; end

    if numel(I)>nSample
        I=I(sort(randperm(numel(I),nSample))) ;
    end

    [du,~]=duvdBGeneral(CtrlVar,MUA,F,l,BCs,KdFuvduv,I) ;
    [~,iMax]=max(vecnorm(du)) ;
    Nodes=[Nodes ; I(iMax)] ;  %#ok<AGROW>

end

end
