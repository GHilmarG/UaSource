function [UserVar,RunInfo,h1,l]=MassContinuityEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%
%
% $$\rho \partial h/\partial t + \partial  (\rho u h)/ \partial x + \partial  (\rho v h)/\partial y - \nabla \cdot (\kappa \nabla h ) = a(h)$$
%
%
%

narginchk(7,8)
nargoutchk(7,7)


[UserVar,RunInfo,h1,l]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);



end



