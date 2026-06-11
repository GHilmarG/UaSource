

function [Kdudh,Kdvdh,Kdhdh,KdFuvhduvh]=duvhdB(CtrlVar,MUA,F,l,BCs,KdFuvhduvh,Nodes)

narginchk(5,7)

if nargin < 6
    KdFuvhduvh=[];
end

if nargin<7 || isempty(Nodes)
    Nodes=1:MUA.Nnodes;
end



%% Calculates the sensitivity matrix duvh/dA
%
%
% The $k$-column contains the response in u and v to a perturbation in $B_k$
%
%
% $$\left[\begin{array}{cccc}
% \partial u_1 /\partial B_1  & \partial u_1 /\partial B_2  & \ldots & \partial u_1 /\partial B_n  \\
% \partial u_2 /\partial B_1  & \partial u_2 /\partial B_2  & \ldots & \partial u_2 /\partial B_n  \\
%              .              &              .              &  .  &    .                          \\
% \partial v_1 /\partial B_1  & \partial v_1 /\partial B_2 & \ldots & \partial v_1 /\partial B_n  \\
% \partial v_2 /\partial B_1  & \partial v_2 /\partial B_2 & \ldots & \partial v_2 /\partial B_n  \\
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
% If, as is generally the case, I have several $q$ variables and several forward models,
%
% $$ F^{uvh}_x =0 $$
%
% $$ F^{uvh}_y =0 $$
%
% $$ F^{uvh}_h =0 $$
%
% then
%
% $$ \frac{\partial F^{uvh}_x}{\partial u} \; \frac{\partial u }{ \partial B}    +  \frac{\partial F^{uvh}_x}{\partial v} \; \frac{\partial v }{ \partial B} + \frac{\partial F^{uvh}_x}{\partial h} \; \frac{\partial h }{ \partial B} = - \frac{\partial F^{uvh}_x }{ \partial B}  $$
%
% $$ \frac{\partial F^{uvh}_y}{\partial u} \; \frac{\partial u }{ \partial B}    +  \frac{\partial F^{uvh}_y}{\partial v} \; \frac{\partial v }{ \partial B} + \frac{\partial F^{uvh}_y}{\partial h} \; \frac{\partial h }{ \partial B} = - \frac{\partial F^{uvh}_y }{ \partial B}  $$
%
% $$ \frac{\partial F^{uvh}_h}{\partial u} \; \frac{\partial u }{ \partial B}    +  \frac{\partial F^{uvh}_h}{\partial v} \; \frac{\partial v }{ \partial B} + \frac{\partial F^{uvh}}{\partial h} \; \frac{\partial h }{ \partial B} = - \frac{\partial F^{uvh}_h }{ \partial B}  $$
%
% or
%
% $$
% \left (\begin{array}{ccc}
% \frac{\partial F^{uvh}_x}{\partial u}  & \frac{F^{uvh}_x}{\partial v}  & \frac{\partial F^{uvh}_x}{\partial h} \\
% \frac{\partial F^{uvh}_y}{\partial u}  & \frac{F^{uvh}_y}{\partial v}  & \frac{\partial F^{uvh}_y}{\partial h} \\
% \frac{\partial F^{uvh}_h} {\partial u}  & \frac{\partial F^{uvh}_h}{\partial v}  & \frac{\partial F^{uvh}_h}{\partial h}
% \end{array}\right )
% \left (\begin{array}{c}
%  \frac{\partial u}{\partial B} \\
%  \frac{\partial v}{\partial B} \\
%  \frac{\partial h}{\partial B}
% \end{array} \right )
% = - \left ( \begin{array}{c}
%   \frac{\partial F^{uvh}_x}{\partial B} \\
%   \frac{\partial F^{uvh}_y}{\partial B} \\
%   \frac{\partial F^{uvh}_h}{\partial B}
%   \end{array} \right )
% $$
%
%
% I assume that where we have boundary conditions on u, v and h, the corresponding sensitives are zero, as any change in p
% (A, B, C) can not affect q (u,v, dh/dt) at those nodes.
%
% This boundary conditions are linear, and the problem itself is linear. So this can be solved as
%
% $$
% \left ( \begin{array}{cc}
%       \frac{\partial F}{\partial q}  & L^T \\
%       L       &  0
% \end{array} \right )
% \left ( \begin{array}{c}
%       \xi  \\
%        0
% \end{array} \right )
% =
% - \left ( \begin{array}{c}
%       \frac{\partial F}{\partial p}  \\
%        0
% \end{array} \right )
% $$
%
%
% see also: duvdAFunc.m, duvdBFunc.m, duvdCFunc.m, dFuvdA.m, dFuvdB.m, dFuvdC.m, TestSensitivityMatrixCalculations.m
%
%%

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;

UserVar=[]; RunInfo=UaRunInfo();

if isempty(KdFuvhduvh)
    [~,~,~,KdFuvhduvh]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F,F,l,BCs) ;
end

%
% KdFuvhduvh=
%    [Kxu Kxv Kxh]
%    [Kyu Kyv Kyh]
%    [Khu Khv Khh]
%
n=MUA.Nnodes;
KdFuvhdh=KdFuvhduvh(1:3*n,2*n+1:3*n);   % Note: for this right-hand-side I'm trying to solve for duvh/dh, not duvh/dB




% BCs
if numel(BCs.ubFixedValue) > 0
    BCs.ubFixedValue=BCs.ubFixedValue*0;
end
if numel(BCs.vbFixedValue) > 0
    BCs.vbFixedValue=BCs.vbFixedValue*0;
end
if numel(BCs.hFixedValue) > 0
    BCs.hFixedValue=BCs.hFixedValue*0;
end

[L,cuvh,luvh]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs);

if ~isempty(L)
    frhs=-KdFuvhdh(:,Nodes) ; 
    grhs=repmat(cuvh,1,size(frhs,2));
else
    frhs=-KdFuvhdh(:,Nodes) ;
    grhs=[];
end

dub=zeros(MUA.Nnodes,1) ; dvb=zeros(MUA.Nnodes,1) ; dh=zeros(MUA.Nnodes,1 ) ; dl=zeros(numel(luvh),1);
sol=solveKApe(KdFuvhduvh,L,frhs,grhs,[dub;dvb;dh],dl,CtrlVar);


Kdudh=sol(1:MUA.Nnodes,:);
Kdvdh=sol(MUA.Nnodes+1:2*MUA.Nnodes,:);
Kdhdh=sol(2*MUA.Nnodes+1:3*MUA.Nnodes,:);



end