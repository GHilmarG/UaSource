

function [dudC,dvdC]=duvdCFunc(CtrlVar,MUA,F,BCs)

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
% $$ \partial F/\partial q \; \partial q / \partial p + \partial F / \partial p = 0 $$
%
% or
%
% $$ \frac{\partial F}{\partial q} \; \frac{\partial q }{ \partial p} = - \frac{\partial F }{ \partial p}  $$
%
% which can be solved for the sensitives
%
% $$ \xi_{ij} : = \frac{\partial q_i}{\partial p_j} $$ 
% 
% If, as is generally the case, I have several $q$ variables then 
%
% $$ \frac{\partial F_u}{\partial u} \; \frac{\partial u }{ \partial C}    +  \frac{\partial F_u}{\partial v} \; \frac{\partial v }{ \partial C} + \frac{\partial F_u}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial C} = - \frac{\partial F_u }{ \partial C}  $$
%
% $$ \frac{\partial F_v}{\partial u} \; \frac{\partial u }{ \partial C}    +  \frac{\partial F_v}{\partial v} \; \frac{\partial v }{ \partial C} + \frac{\partial F_v}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial C} = - \frac{\partial F_v }{ \partial C}  $$
%
% $$ \frac{\partial F_{\dot{h}}}{\partial u} \; \frac{\partial u }{ \partial C}    +  \frac{\partial F_{\dot{h}}}{\partial v} \; \frac{\partial v }{ \partial C} + \frac{\partial F_{\dot{h}}}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial C} = - \frac{\partial F_{\dot{h}} }{ \partial C}  $$
%
%
% Note: It is here assumed that the forward problem has already been solved. So ahead of a call to this function one needs to
% have called
% 
%  [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
%
%
% and the F provided as an input to this function must be this solution to the forward problem.
%
% see also: dFuvdA.m, dFuvdC.m , TestSensitivityMatrixCalculations.m
% 
%%

dFdC=dFuvdC(CtrlVar,MUA,F) ;
dFdC=-dFdC; % there is actually a different sign convention inside of this...

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;

[~,dFduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);

% test for one particular node and compare to finite-differences 

% if velocities are prescribed, the sensitivity of those velocities to changes in model parameters is zero.
% make sure that the BCs reflect this.
if numel(BCs.ubFixedValue) > 0
    BCs.ubFixedValue=BCs.ubFixedValue*0;
end
if numel(BCs.vbFixedValue) > 0
    BCs.vbFixedValue=BCs.vbFixedValue*0;
end

[L,cuv]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;

if isempty(cuv)
    l.ubvb=[];
else
    l.ubvb=zeros(numel(cuv),1) ;
end

if ~isempty(L)
    frhs=-dFdC-L'*l.ubvb; % Note, this uses Matlab automatic implicit expansion to expand the L'*l column to match the dimensions of the dFdA matrix
    grhs=cuv-L*[F.ub;F.vb] ;
else
    frhs=-dFdC ;
    grhs=[];
end

duvb=zeros(2*MUA.Nnodes,1) ; dl=zeros(numel(l.ubvb),1);
CtrlVar.TestKApeSolve=false; 
sol=solveKApe(dFduv,L,frhs,grhs,duvb,dl,CtrlVar);

%l.ubvb=dl; 

dudC=sol(1:MUA.Nnodes,:);
dvdC=sol(MUA.Nnodes+1:end,:);



end