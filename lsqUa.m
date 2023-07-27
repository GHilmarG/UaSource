function [x,lambda,R2,Slope0,g2,residual,g,h,output] = lsqUa(CtrlVar,fun,x,lambda,L,c)

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

if isempty(CtrlVar) || ~isstruct(CtrlVar)

    Algorithm="DogLeg";

else

    Algorithm=CtrlVar.lsqUa.Algorithm;

end



switch Algorithm

    case "LevenbergMarquardt"

        [x,lambda,R2,Slope0,g2,residual,g,h,output] = lsqLevenbergMarquardtUa(CtrlVar,fun,x,lambda,L,c) ;


    case "DogLeg"


        [x,lambda,R2,Slope0,g2,residual,g,h,output] = lsqDogLegUa(CtrlVar,fun,x,lambda,L,c) ;


    otherwise

        error("case not found")

end