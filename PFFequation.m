function [UserVar,phi,lambda,HEmatrix,HErhs]=PFFequation(UserVar,CtrlVar,MUA,BCs,Gc,l,Psi)

%% 
%
% Solves
%
% $$  G_c \left ( \frac{\phi(x,y)}{l} - l \nabla^2 \phi(x,y) \right ) = (1-\phi) \Psi(x,y) $$
%
% By solving the in-homogeneous Helmholtz equation with variable coefficients in two dimensions:
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y) - \nabla \cdot \nabla d(x,y) $$
%
% or
%
% $$  (G_c/l +\Psi(x,y) ) \, \phi(x,y) - G_c \, l \, \nabla \cdot \nabla \phi(x,y)) = \Psi(x,y) $$
%
%
% with
%
% $$a(x,y)=G_c/l+\Psi(x,y) $$
%
% $$b=G_c \, l$$
%
% $$ c(x,y)=\Psi(x,y) $$
%
% $$ d=0 $$
%
%%
narginchk(7,7)

N=MUA.Nnodes; 
a=Gc/l+Psi+zeros(N,1);
b=Gc*l+zeros(N,1);
c=Psi+zeros(N,1) ;
d=zeros(N,1);

[UserVar,HEmatrix,HErhs]=HelmholtzEquationAssembly(UserVar,CtrlVar,MUA,a,b,c,d);



% Now apply BCs.  
% Note: When defining tracer boundary conditions, use the thickness fields (h) in the BCs structure
% for that purpose.
%
MLC=BCs2MLC(CtrlVar,MUA,BCs);
L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;


phi0=zeros(N,0) ;
[phi,lambda]=solveKApe(HEmatrix,L,HErhs,Lrhs,phi0,lambda,CtrlVar);
phi=full(phi);






end




