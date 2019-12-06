function [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d,RHS)

%%
% Solves the Helmholtz equation:
%
% $$  \kappa^2 f(x,y) - \nabla^2 f(x,y)) = c(x,y) - \nabla d(x,y) $$
%
%
% The non-homogeneous Helmholtz equaton with variable coefficents in two spatial
% dimentions is
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y) - \nabla \cdot \nabla  d(x,y) $$
%
% which we can also write as
%
% $$  a ( f - \tilde{f} ) - \nabla \cdot (b \nabla (f-\tilde{f})) = 0 $$
%
% where $\tilde{f}$ is a given function

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




