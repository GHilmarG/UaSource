function [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,RHS)

%%
% Solves the Helmholtz equation:
%
% $$  \kappa^2 f(x,y) + \nabla^2 f(x,y)) = a$$
%
%
% The non-homogeneous Helmholtz equaton with variable coefficents in two spatial
% dimentions is
%
% $$  a(x,y) f(x,y) - \nabla \cdot (b(x,y) \nabla f(x,y)) = c(x,y)$$
%

[UserVar,HEmatrix,HErhs]=HelmholtzEquationAssembly(UserVar,CtrlVar,MUA,a,b,c);
L=[] ; Lrhs=[] ; lambda=[]; f=[] ;

if nargin==7 && ~isempty(RHS)
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




