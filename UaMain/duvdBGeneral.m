


function [dudB,dvdB]=duvdBGeneral(CtrlVar,MUA,F,l,BCs,KdFuvduv,Nodes)  %#ok<INUSL>

%% Calculates the sensitivity matrices du/dB and dv/dB
%
% If $n$ is the number of nodes, the matrices returned are each $n \times n$, and the $k$-th column contains the
% response in u (or v) to a perturbation in $B_k$
%
% $$\left[\begin{array}{cccc}
% \partial u_1 /\partial B_1  & \partial u_1 /\partial B_2  & \ldots & \partial u_1 /\partial B_n  \\
% \partial u_2 /\partial B_1  & \partial u_2 /\partial B_2  & \ldots & \partial u_2 /\partial B_n  \\
%              .              &              .              &  .     &    .                        \\
% \end{array}\right] $$
%
% This is the general version of
%
%   duvdBFunc.m
%
% Unlike duvdBFunc.m, no assumption is made about the ice being grounded: the floating mask is itself a function of B
% and its variation with B is accounted for, through dFuvdBGeneral.m
%
% Note that, unlike duvdBFunc.m, this function does NOT return the sensitivities of $\dot{h}$ to B. Only the
% sensitivities of the two velocity components are calculated.
%
%% Approach
%
% If the forward model is
%
% $$ F(q(p),p) = 0 $$
%
% where $q$ are the output variables and $p$ the model parameters, then
%
% $$ \frac{\partial F}{\partial q} \; \frac{\partial q }{ \partial p} = - \frac{\partial F }{ \partial p}  $$
%
% which is solved here for the sensitivities, with $q=(u,v)$ and $p=B$.
%
% $\partial F / \partial q$ is the Newton (tangent) matrix of the uv system, and is either provided as an input
% (KdFuvduv) or calculated here. $\partial F / \partial p$ is provided by dFuvdBGeneral.m
%
%% Assumptions
%
% As in dFuvdBGeneral.m and dIdBqGeneral.m: the upper ice surface s, the ocean surface S and the densities are held
% fixed, and b and h are obtained from s, S and B using the closure solved by Calc_bh_From_sBS.m. F.b and F.h must be
% the converged output of that closure for the current F.B.
%
% It is here assumed that the forward problem has already been solved, i.e. that ahead of a call to this function one
% has called
%
%  [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
%
% and that the F provided as input is that solution.
%
%% Inputs
%
%   Nodes    optional. If provided, only the columns of the sensitivity matrices corresponding to those nodes are
%            calculated, i.e. only perturbations in B at those nodes are considered. The returned matrices are then
%            n x numel(Nodes). This can give a very significant saving, both in time and in memory, because the
%            right-hand side of the linear system is dense.
%
%  see also: dFuvdBGeneral.m, dIdBqGeneral.m, dGeometrydB.m, duvdBFunc.m, duvdCFunc.m
%
%%

narginchk(5,7)
nargoutchk(2,2)

if nargin<7 || isempty(Nodes)
    Nodes=1:MUA.Nnodes;
end

if nargin<6 || isempty(KdFuvduv)
    CtrlVar.uvAssembly.ZeroFields=false;
    CtrlVar.uvMatrixAssembly.Ronly=false;
    [~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
end

%%  dF/dp

KdFuvdB=dFuvdBGeneral(CtrlVar,MUA,F);

KdFuvdB=KdFuvdB(:,Nodes);    % only the requested columns are needed

%% boundary conditions
%
% Where velocities are prescribed, the sensitivity of those velocities to changes in the model parameters is zero, so
% the boundary values must be set to zero here.

if numel(BCs.ubFixedValue) > 0
    BCs.ubFixedValue=BCs.ubFixedValue*0;
end

if numel(BCs.vbFixedValue) > 0
    BCs.vbFixedValue=BCs.vbFixedValue*0;
end

[LBCs,cBCs]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;

%% solve

frhs=-full(KdFuvdB);   % the right-hand side is quite dense, so this is a faster approach

if ~isempty(LBCs)
    grhs=repmat(cBCs,1,size(frhs,2));
else
    grhs=[];
end

CtrlVar.TestKApeSolve=false;
sol=solveKApe(KdFuvduv,LBCs,frhs,grhs,[],[],CtrlVar);

dudB=sol(1:MUA.Nnodes,:);
dvdB=sol(MUA.Nnodes+1:2*MUA.Nnodes,:);

end
