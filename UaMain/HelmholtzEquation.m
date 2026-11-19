

function [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d,RHS,BCs)

%% 
% Solves the in-homogeneous Helmholtz equation with variable coefficients in two dimensions:
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y) - \nabla \cdot \nabla d(x,y) $$
%
% Also possible to specify the right-hand-side directly through 
%
%   RHS
%
% If RHS is given as input, c and d are not used. 
%
% The in-homogeneous Helmholtz equation with variable coefficients in two spatial
% dimensions is
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y) - \nabla \cdot \nabla  d(x,y) $$
%
% which we can also write as
%
% $$  a ( f - \tilde{f} ) - \nabla \cdot (b \nabla (f-\tilde{f})) = 0 $$
%
% where  $c= -a \tilde{f}$ and $d= -b \tilde{f}$, and $\tilde{f}$ is a given function
%
% Boundary conditions:
%
% Without BCs (or with BCs=[]) natural boundary conditions apply (as before). Dirichlet conditions, and ties between nodes, can be
% given in BCs, a BoundaryConditions object, of which the thickness fields are used:
%
%   BCs.hFixedNode, BCs.hFixedValue   : nodes where f is prescribed, and the values
%   BCs.hTiedNodeA, BCs.hTiedNodeB    : f(hTiedNodeA)=f(hTiedNodeB) (eg periodic boundary conditions)
%   (BCs.hPosNode, BCs.hPosValue are treated as BCs.hFixedNode, BCs.hFixedValue)
%
% The constraints are formed by BCs2MLC.m (MLC.hL, MLC.hRhs), as elsewhere in Ua, and imposed with Lagrange multipliers
% (solveKApe.m), which are returned in lambda. BCs2MLC.m checks that no node is constrained more than once.
%
% Examples:  
%
% Smooth a given field over a FE mesh:
%
%    load('PIG-TWG-RestartFile.mat') ; CtrlVar=CtrlVarInRestartFile;
%    L=1e3 ;  % Smoothing length scale 
%    [UserVar,SmoothedField]=HelmholtzEquation([],CtrlVar,MUA,1,L^2,F.B,0); 
%
%    figure(1) ; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,SmoothedField) ; title(' Smoothed field') ; xlabel('x (km)') ; ylabel('y (km)') 
%    figure(2) ; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,F.B) ; title(' Original field')  ; xlabel('x (km)') ; ylabel('y (km)') 
%    figure(3) ; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,SmoothedField-F.B) ; title(' Smoothed-Original')  ; xlabel('x (km)') ; ylabel('y (km)') 
%
% A field that is prescribed at the boundary nodes and decays inwards over the length L (Dirichlet conditions):
%
%    BCs=BoundaryConditions; BCs.hFixedNode=MUA.Boundary.Nodes; BCs.hFixedValue=ValuesAtTheBoundary;
%    [UserVar,f]=HelmholtzEquation([],CtrlVar,MUA,1,L^2,0,0,[],BCs);
%
%
%%
narginchk(7,9)


[UserVar,HEmatrix,HErhs]=HelmholtzEquationAssembly(UserVar,CtrlVar,MUA,a,b,c,d);
L=[] ; Lrhs=[] ; lambda=[]; f=[] ;

if nargin>=8 && ~isempty(RHS)
    HErhs=RHS;
end

% Dirichlet conditions and ties (the h fields of BCs), imposed with Lagrange multipliers, through BCs2MLC.m as elsewhere in Ua
if nargin==9 && ~isempty(BCs)
    MLC=BCs2MLC(CtrlVar,MUA,BCs);
    L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;
end

[f,lambda]=solveKApe(HEmatrix,L,HErhs,Lrhs,f,lambda,CtrlVar);
f=full(f);


end
