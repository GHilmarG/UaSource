function [UserVar,c1,lambda]=LevelSetEquation(UserVar,CtrlVar,MUA,F,phi0,SF)

%%
% Solves the tracer conservation equation for the tracer c on the form:
% 
% $$\partial c/\partial t + d (u c)/dx + d (v c)/dy - \nabla \cdot (\kappa \nabla c) = a$$
% 
% The natural boundary condition is 
%
% $$\nabla c \cdot \hat{u} = 0 $$
%
% ie,  the free outflow condition 
%


%   < f | N + M >  

dx=phi0*Dxx*phiy
dy=;phi0*Dyy*phiy
Normphie=vecnorm([dx(:) dy(:)],2,2)


MLC=BCs2MLC(CtrlVar,MUA,BCsTracer);
L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;

[UserVar,kv,rh]=TracerConservationEquationAssembly(UserVar,CtrlVar,MUA,dt,c0,u0,v0,a0,u1,v1,a1,kappa);

[c1,lambda]=solveKApe(kv,L,rh,Lrhs,c0,lambda,CtrlVar);
c1=full(c1);


end

