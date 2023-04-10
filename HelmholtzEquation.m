function [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d,RHS)

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
% The in-homogeneous Helmholtz equation with variable coefficents in two spatial
% dimentions is
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y) - \nabla \cdot \nabla  d(x,y) $$
%
% which we can also write as
%
% $$  a ( f - \tilde{f} ) - \nabla \cdot (b \nabla (f-\tilde{f})) = 0 $$
%
% where  $c= -a \tilde{f}$ and $d= -b \tilde{f}$, and $\tilde{f}$ is a given function
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
%
%%
narginchk(7,8)


[UserVar,HEmatrix,HErhs]=HelmholtzEquationAssembly(UserVar,CtrlVar,MUA,a,b,c,d);
L=[] ; Lrhs=[] ; lambda=[]; f=[] ;

if nargin==8 && ~isempty(RHS)
    HErhs=RHS;
end

[f,lambda]=solveKApe(HEmatrix,L,HErhs,Lrhs,f,lambda,CtrlVar);
f=full(f);

%
% MLC=BCs2MLC(CtrlVar,MUA,BCsTracer);
% L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;
% [c1,lambda]=solveKApe(kv,L,rh,Lrhs,c0,lambda,CtrlVar);
% c1=full(c1);


end




