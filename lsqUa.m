



function [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x,lambda,L,c)

%%
%
% Minimizes the norm of R where R is a vector subject the the
% constraints
%
%   L x = c
%
%  The function 'fun' should provide both R and the Jacobian, i.e.
%
%   [R,K]=fun(x)
%
% where   K=grad J
%
%
%  fun = @(x) FunctionThatReturnsRandK(x, ...other variables whos value do not change as the function is called repeatedly...)  ;
%
%
% Note that R is a vector and K a matrix.  The value that is minimized is
%
%   R'*R
%
% Algorithm: Solves repeatedly the linearized min problem:
%
%  |R0+J dx|^2
%
%
% Example:
%
%   lsqUaExample
%
%%

if isempty(CtrlVar) || ~isstruct(CtrlVar) || ~isfield(CtrlVar,"lsqUa") || ~isfield(CtrlVar.lsqUa,"Algorithm")

    Algorithm="DogLeg";

else

    Algorithm=CtrlVar.lsqUa.Algorithm;

end



switch Algorithm

    case "LevenbergMarquardt"

        [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqLevenbergMarquardtUa(CtrlVar,fun,x,lambda,L,c) ;
        

    case "DogLeg"


        [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqDogLegUa(CtrlVar,fun,x,lambda,L,c) ;


    otherwise

        error("case not found")

end