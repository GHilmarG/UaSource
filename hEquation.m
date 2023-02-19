function [UserVar,h,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs,kIso,kAlong,kCross)

%%
% Solves the linear h-equation on the form:
% 
% $$ d (u h)/dx + d (v h)/dy - \nabla \cdot (\kappa \nabla h) = a$$
% 
% for h
%
% The natural boundary condition is 
%
%  $$\nabla h \cdot \hat{n} = 0 $$
%
% ie,  the free outflow condition 
%
% Boundary conditions: The BCs are identical to define thickness boundary
% conditions. 
%


[UserVar,Kh,bh]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as,kIso,kAlong,kCross);

% Now apply BCs.  
% Note: When defining tracer boundary conditions, use the thickness fields (h) in the BCs structure
% for that purpose.
%
MLC=BCs2MLC(CtrlVar,MUA,BCs);
L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;


h0=F.x*0;
[h,lambda]=solveKApe(Kh,L,bh,Lrhs,h0,lambda,CtrlVar);
h=full(h);


end




