

function [dudC,dvdC,dhdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv,Nodes)

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
narginchk(5,7)

if nargin<7 || isempty(Nodes)
    Nodes=1:MUA.Nnodes;
end

if nargin < 6 || isempty(KdFuvduv)
    CtrlVar.uvAssembly.ZeroFields=false;
    CtrlVar.uvMatrixAssembly.Ronly=false;
    [~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);
end

Sensitivities="-dudC-dvdC-";

% if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
%     Sensitivities=Sensitivities+"-dhdotdC-" ;
%     Sensitivities=replace(Sensitivities,"--","-");
% end
% 

KdFuvdC=dFuvdC(CtrlVar,MUA,F) ;


if numel(BCs.ubFixedValue) > 0
    BCs.ubFixedValue=BCs.ubFixedValue*0;
end
if numel(BCs.vbFixedValue) > 0
    BCs.vbFixedValue=BCs.vbFixedValue*0;
end


if contains(Sensitivities,"dhdotdC")
    if numel(BCs.hFixedValue) > 0
        BCs.hFixedValue=BCs.hFixedValue*0;
    end
    [LBCs,cBCs]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs);
    [KdFhdotduvhdot]=dFhdot_duvhdot(CtrlVar,MUA,F) ;
    O2n1n=sparse(2*MUA.Nnodes,MUA.Nnodes);
    KdFdq=[KdFuvduv O2n1n ; KdFhdotduvhdot];
    Onn=sparse(MUA.Nnodes,MUA.Nnodes);
    KdFdp=[KdFuvdC ; Onn] ;

else
    KdFdq=KdFuvduv ; 
    KdFdp=KdFuvdC;
    [LBCs,cBCs]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;
end


if ~isempty(LBCs)
    frhs=-KdFdp ;
    grhs=repmat(cBCs,1,size(frhs,2));
else
    frhs=-KdFdp ;
    grhs=[];
end

CtrlVar.TestKApeSolve=false;
sol=solveKApe(KdFdq,LBCs,frhs,grhs,[],[],CtrlVar);


dudC=sol(1:MUA.Nnodes,:);
dvdC=sol(MUA.Nnodes+1:2*MUA.Nnodes,:);
if contains(Sensitivities,"dhdotdC")
    dhdC=sol(2*MUA.Nnodes+1:3*MUA.Nnodes,:);
else
    dhdC=[];
end



end