function [UserVar,RunInfo,h1,l]=MassContinuityEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%  $$\rho \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \mathbf{v} h )  =  a(h)$$ 
%
%
%%

narginchk(8,8)
nargoutchk(4,4)


[UserVar,RunInfo,h1,l]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);



end



