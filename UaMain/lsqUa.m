



function [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x,lambda,L,c)

%%
%
% Minimizes the norm or R where R is a vector subject to the
% constraints
%
%   $$L x = c$$
%
%  The function 'fun' should provide both R and the Jacobian, i.e.
%
%   [R,K]=fun(x)
%
% where   K=grad R
%
%
%  fun = @(x) FunctionThatReturnsRandK(x, ...other variables which values do not change as the function is called repeatedly...)  ;
%
%
% Note that R is a vector and K a matrix.  The value that is minimized is
%
%   R'*R
%
% *Algorithm:* 
% 
% Solves repeatedly the linearized min problem:
%
%  |R0+J dx|^2
%
%
%
% The KKT systems are:
%
% Newton:
%
%   [K0 L' ]  [dx]  = - [R0 - L' lambda0]
%   [L  0  ]  [dl]      [L x0 - c        ]
%
%
% Newton lsq
%
%   [2 K0'K0   L' ]  [dx]  = - [2 K0' R0 + L' lambda0 ]
%   [    L     0  ]  [dl]      [      L x0 - c        ]
%
%
% Cauchy:
%
%   [    I     L' ]  [dx]  = - [2 K0' R0 + L' lambda0 ]
%   [    L     0  ]  [dl]      [      L x0 - c        ]
%
% search direction:  d=[dx ; dlambda]
%
% Algorithm: Solves a least-squares problem such as
%
% 
% $$\min_{x}  R^2 = \| \mathbf{R} \|^2 $$
% 
% with $\mathbf{R}$ being a vector.
%
% We assume we know
%
% $$K=\nabla \mathbf{R}$$
%
%
% And we find that 
%
% $$ \nabla R = K'\mathbf{R}  $$
%
% and
%
% $$ \nabla^2 R  \approx K' K$$
%
%
% The quadratic approximation is therefore
%
% $$ R \approx Q := R_0 + K' \mathbf{R} \, \Delta x + \frac{1}{2} \Delta x \, K' K \Delta x $$
%
%
% and the Newton system is 
%
% $$ K' K \Delta x =  -K' \mathbf{R}$$  
%
% If $K$ is $n \times n$ and invertible, this is same as solving
%
% $$ K \Delta x = - \mathbf{R} $$
%
% and gives the same update and direction.
%
% So the Newton system is the same. 
% 
% However, when finding the Cauchy step we must search in the direction
%
% $$ K' \mathbf{R} $$
%
% and not simply along $\mathbf{R}$ 
%
% Various exit criteria are possible.
%
% Parameters: 
%
%   CtrlVar.lsqUa.Algorithm=["LevenbergMarquardt"|"DogLeg"]
%   CtrlVar.lsqUa.ItMax=5; 
%   CtrlVar.lsqUa.gTol=1e-20;
%   CtrlVar.lsqUa.dR2Tol=1e-3; 
%   CtrlVar.lsqUa.dxTol=1e-20;
%   CtrlVar.lsqUa.isLSQ=true; 
%   CtrlVar.lsqUa.SaveIterate=false; 
%   CtrlVar.lsqUa.CostMeasure=["r2"|"R2"]  
%
% r2 is the norm of the right-hand side of the KKT system.
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