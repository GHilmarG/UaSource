

function [dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs,Nodes)

%% Calculates the sensitivity matrix duv/dA
%
% If $n$ is the number of nodes, the matrix returned is $2 n \times n$ 
%
% The $k$-column contains the response in u and v to a perturbation in $A_k$
%
%
% $$\left[\begin{array}{cccc} 
% \partial u_1 /\partial A_1  & \partial u_1 /\partial A_2  & \ldots & \partial u_1 /\partial A_n  \\  
% \partial u_2 /\partial A_1  & \partial u_2 /\partial A_2  & \ldots & \partial u_2 /\partial A_n  \\  
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
% $$ \frac{\partial F_u}{\partial u} \; \frac{\partial u }{ \partial A}    +  \frac{\partial F_u}{\partial v} \; \frac{\partial v }{ \partial A} + \frac{\partial F_u}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial A} = - \frac{\partial F_u }{ \partial A}  $$
%
% $$ \frac{\partial F_v}{\partial u} \; \frac{\partial u }{ \partial A}    +  \frac{\partial F_v}{\partial v} \; \frac{\partial v }{ \partial A} + \frac{\partial F_v}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial A} = - \frac{\partial F_v }{ \partial A}  $$
%
% $$ \frac{\partial F_{\dot{h}}}{\partial u} \; \frac{\partial u }{ \partial A}    +  \frac{\partial F_{\dot{h}}}{\partial v} \; \frac{\partial v }{ \partial A} + \frac{\partial F_{\dot{h}}}{\partial \dot{h}} \; \frac{\partial \dot{h} }{ \partial A} = - \frac{\partial F_{\dot{h}} }{ \partial A}  $$
%
% or
%
% $$
% \left (\begin{array}{ccc} 
% \frac{F_u}{\partial u}  & \frac{F_u}{\partial v}  & \frac{\partial F_u}{\partial \dot{h}} \\
% \frac{F_v}{\partial u}  & \frac{F_v}{\partial v}  & \frac{\partial F_v}{\partial \dot{h}} \\
% \frac{F_{\dot{h}}} {\partial u}  & \frac{F_{\dot{h}}}{\partial v}  & \frac{\partial F_{\dot{h}}}{\partial \dot{h}} 
% \end{array}\right ) 
% \left (\begin{array}{c} 
%  \frac{\partial u}{\partial A} \\
%  \frac{\partial v}{\partial A} \\ 
%  \frac{\partial \dot{h}}{\partial A} 
% \end{array} \right )
% = - \left ( \begin{array}{c}
%   \frac{\partial F_u}{\partial A} \\
%   \frac{\partial F_v}{\partial A} \\
%   \frac{\partial F_{\dot{h}}}{\partial A}  
%   \end{array} \right )
% $$
%
% The momentum equations are not explicit functions of $\dot{h}$ and the mass-conservation equation is not a function of $A$.
% Therefore the system above simplifies to
%
% $$
% \left (\begin{array}{ccc} 
% \frac{F_u}{\partial u}  & \frac{F_u}{\partial v}  & 0 \\
% \frac{F_v}{\partial u}  & \frac{F_v}{\partial v}  & 0 \\
% \frac{F_{\dot{h}}} {\partial u}  & \frac{F_{\dot{h}}}{\partial v}  & I
% \end{array}\right ) 
% \left (\begin{array}{c} 
%  \frac{\partial u}{\partial A} \\
%  \frac{\partial v}{\partial A} \\ 
%  \frac{\partial \dot{h}}{\partial A} 
% \end{array} \right )
% = - \left ( \begin{array}{c}
%   \frac{\partial F_u}{\partial A} \\
%   \frac{\partial F_v}{\partial A} \\
%   0
%   \end{array} \right )
% $$
%
% Note that $I$ is here not the identity matrix but the mass matrix.
%
%
% One could therefore also first solve 
%
% $$
% \left (\begin{array}{cc} 
% \frac{F_u}{\partial u}  & \frac{F_u}{\partial v} \\
% \frac{F_v}{\partial u}  & \frac{F_v}{\partial v}  
% \end{array}\right ) 
% \left (\begin{array}{c} 
%  \frac{\partial u}{\partial A} \\
%  \frac{\partial v}{\partial A} 
% \end{array} \right )
% = - \left ( \begin{array}{c}
%   \frac{\partial F_u}{\partial A} \\
%   \frac{\partial F_v}{\partial A} 
%   \end{array} \right )
% $$
%
% and then determine $\partial{\dot{h}}/\partial{A}$ from
%
% $$
%  \frac{\partial \dot{h}}{\partial A} = 
% - \frac{\partial F_{\dot{h}}}{\partial u} \frac{\partial u}{\partial A} 
% - \frac{\partial F_{\dot{h}}}{\partial v} \frac{\partial v}{\partial A} 
% $$
%
% Note: It is here assumed that the forward problem has already been solved. So ahead of a call to this function one needs to
% have called
% 
%  [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
%
%
% and the F provided as an input to this function must be this solution to the forward problem.
%
% see also: duvdAFunc.m, duvdBFunc.m, duvdCFunc.m, dFuvdA.m, dFuvdB.m, dFuvdC.m, TestSensitivityMatrixCalculations.m
% 
%%

if nargin<5 || isempty(Nodes)
    Nodes=1:MUA.Nnodes;
end


KdFuvdA=dFuvdA(CtrlVar,MUA,F);


CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;

[~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


%%  To do for dh/dt meas
%
%   [KdFhdotdu,KdFhdotdv,KdFhdotdhdot]=dFhdot_duvhdot(CtrlVar,MUA,F) ; 
%
% then solve:
%
%   [       KdFuvdu           0          ] [duv/dA ]  = - [ KdFuvdA   ]
%   [KdFhdotdu KdFhdotdv KdFhdothdot     ] [dhdot/dA]     [ KdFhdotdA ] 
%
%
% where 
%
%  KdFhdothdot = M 
%  KdFhdotA    = 0 
%
% and need to add in BCs for u,v, and h
%
%%


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
    frhs=-KdFuvdA(:,Nodes)-L'*l.ubvb; % Note, this uses Matlab automatic implicit expansion to expand the L'*l column to match the dimensions of the dFdA matrix
    %frhs=-dFdA(:,Node)-L'*l.ubvb; % if only calculate for one given node
    grhs=cuv-L*[F.ub;F.vb] ;
else
    frhs=-KdFuvdA ;
    grhs=[];
end


dub=zeros(MUA.Nnodes,1) ; dvb=zeros(MUA.Nnodes,1) ; dl=zeros(numel(l.ubvb),1);
CtrlVar.TestKApeSolve=false; 
sol=solveKApe(KdFuvduv,L,frhs,grhs,[dub;dvb],dl,CtrlVar);


dudA=sol(1:MUA.Nnodes,:);
dvdA=sol(MUA.Nnodes+1:end,:);


end