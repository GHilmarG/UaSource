

function duvC=duvdCFunc(CtrlVar,MUA,F,BCs)

%% Calculates the sensitivity matrix duv/dC
%
% If $n$ is the number of nodes, the matrix returned is $2 n \times n$ 
%
% The $k$-column contains the response in u and v to a perturbation in $C_k$
%
%
% $$\left[\begin{array}{cccc} 
% \partial u_1 /\partial C_1  & \partial u_1 /\partial C_2  & \ldots & \partial u_1 ;\partial C_n  \\  
% \partial u_2 /\partial C_1  & \partial u_2 /\partial C_2  & \ldots & \partial u_2 ;\partial C_n  \\  
%              .              &              .              &  .  &    .                          \\  
% \partial v_1 /\partial C_1  & \partial v_1 /\partial C_2 & \ldots & \partial v_1 /\partial C_n  \\
% \partial v_2 /\partial C_1  & \partial v_2 /\partial C_2 & \ldots & \partial v_2 /\partial C_n  \\
%              .              &              .              &  .  &    .                          \\  
% \end{array}\right] $$
%
% Approach:
%
% If the forward model is
%
% $$ F(q(p),p) = 0 $$
%
% where $q$ are output variables and $p$ model parameters, then
%
% $$ \partial F/\partial q \; \partial q / \partial q + \partial F / \partial p = 0 $$
%
% or
%
% $$ \partial F/\partial q \; \partial q / \partial q = - \partial F / \partial p  $$
%
% which can be solved for the sensitives
%
% $$ \partial q / \partial q $$ 
% 
% Note: It is here assumed that the forward problem has already been solved. So ahead of a call to this function one needs to
% have called
% 
%  [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
%
%
% and the F provided as an input to this function must be this solution to the forward problem.
%
%%


dFdp=dFuvdC(CtrlVar,MUA,F) ;

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;

[~,dFdq]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


duvC=dFdq\dFdp; % Note a minus is not missing because of the way dFuvdC is constructed (might want to change this later...)

end