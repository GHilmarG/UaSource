

function [dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs)

%% Calculates the sensitivity matrix duv/dA
%
% If $n$ is the number of nodes, the matrix returned is $2 n \times n$ 
%
% The $k$-column contains the response in u and v to a perturbation in $A_k$
%
%
% $$\left[\begin{array}{cccc} 
% \partial u_1 /\partial A_1  & \partial u_1 /\partial A_2  & \ldots & \partial u_1 ;\partial A_n  \\  
% \partial u_2 /\partial A_1  & \partial u_2 /\partial A_2  & \ldots & \partial u_2 ;\partial A_n  \\  
%              .              &              .              &  .  &    .                          \\  
% \partial v_1 /\partial A_1  & \partial v_1 /\partial A_2 & \ldots & \partial v_1 /\partial A_n  \\
% \partial v_2 /\partial A_1  & \partial v_2 /\partial A_2 & \ldots & \partial v_2 /\partial A_n  \\
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
% see also: duvdCFunc.m, dFuvdA.m, dFuvdC.m , TestSensitivityMatrixCalculations.m
% 
%%

dFdA=dFuvdA(CtrlVar,MUA,F);

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;

[~,dFduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


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
    frhs=-dFdA-L'*l.ubvb; % Note, this uses Matlab automatic implicit expansion to expand the L'*l column to match the dimensions of the dFdA matrix
    %frhs=-dFdA(:,Node)-L'*l.ubvb; % if only calculate for one given node
    grhs=cuv-L*[F.ub;F.vb] ;
else
    frhs=-dFdA ;
    grhs=[];
end


dub=zeros(MUA.Nnodes,1) ; dvb=zeros(MUA.Nnodes,1) ; dl=zeros(numel(l.ubvb),1);
CtrlVar.TestKApeSolve=false; 
sol=solveKApe(dFduv,L,frhs,grhs,[dub;dvb],dl,CtrlVar);


dudA=sol(1:MUA.Nnodes,:);
dvdA=sol(MUA.Nnodes+1:end,:);


end