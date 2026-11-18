

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
% Boundary conditions (18 Nov 2026):
%
% Without BCs (or with BCs=[]) natural boundary conditions apply (as before). Dirichlet conditions, and ties between nodes, can be
% given in BCs, a BoundaryConditions object or a structure, of which the fields
%
%   BCs.hFixedNode, BCs.hFixedValue   : nodes where f is prescribed, and the values
%   BCs.hTiedNodeA, BCs.hTiedNodeB    : f(hTiedNodeA)=f(hTiedNodeB) (eg periodic boundary conditions)
%
% are used (missing fields are treated as empty). The constraints are imposed with Lagrange multipliers, as elsewhere in Ua
% (createLh.m, ScaleL.m, solveKApe.m), and the multipliers are returned in lambda. A node must not be constrained more than once.
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

% (18 Nov 2026) Dirichlet conditions and ties, imposed with Lagrange multipliers
if nargin==9 && ~isempty(BCs)
    FixedNode=BCsField(BCs,"hFixedNode"); FixedValue=BCsField(BCs,"hFixedValue");
    TiedA=BCsField(BCs,"hTiedNodeA"); TiedB=BCsField(BCs,"hTiedNodeB");
    Constrained=[FixedNode(:);TiedA(:);TiedB(:)];
    if numel(unique(Constrained))<numel(Constrained)
        error("HelmholtzEquation:DuplicateConstraints","A node must not be constrained more than once (in hFixedNode, hTiedNodeA and hTiedNodeB).")
    end
    [L,Lrhs]=createLh(MUA.Nnodes,FixedNode(:),FixedValue(:),TiedA(:),TiedB(:));
    if ~isempty(L)
        [L,Lrhs]=ScaleL(CtrlVar,L,Lrhs);
        lambda=zeros(size(L,1),1);
    end
end

[f,lambda]=solveKApe(HEmatrix,L,HErhs,Lrhs,f,lambda,CtrlVar);
f=full(f);

%
% MLC=BCs2MLC(CtrlVar,MUA,BCsTracer);
% L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;
% [c1,lambda]=solveKApe(kv,L,rh,Lrhs,c0,lambda,CtrlVar);
% c1=full(c1);


end


function v=BCsField(BCs,Name)
% field or property of BCs, empty if missing
if (isstruct(BCs) && isfield(BCs,Name)) || (isobject(BCs) && isprop(BCs,Name))
    v=BCs.(Name) ;
else
    v=[] ;
end
end
