


function H = CalcDirectAdjointHessian(UserVar,CtrlVar,RunInfo,MUA,F,BCs,l,BCsAdjoint,d2Iduu,d2Idvv,d2Idhdothdot,uAdjoint,vAdjoint)

narginchk(13,13)


%% Calculates some of the terms of the Hessian, H, using the direct-adjoint approach.
%
% Here this is done for the misfit term and later the Hessian of the regularization term is added (this is easy).
%
% The tricky part is to do this for the misfit term.
%
%
% $$
%   H_{ij} = \frac{\partial^2 J}{\partial p_i \, \partial p_j}
%   + \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial p_j}
%   +\frac{\partial^2 J}{\partial q_k\, \partial q_m} \xi_{ki} \, \xi_{mj}
%   +\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial q_m} \xi_{ki} \, \xi_{mj}
%   +\frac{\partial^2 J}{\partial p_i \, \partial q_k} \xi_{kj}
%   +\Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj}
%   +\frac{\partial^2 J}{\partial q_k \, \partial p_j} \xi_{ki}
%   +\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki}
% $$
%
%
% The term
%
% $$
%   \frac{\partial^2 J}{\partial p_i \, \partial p_j}
% $$
%
% is easy, and is done in Regulaisation.m
%
%
% The term
%
% $$
%   \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial p_j}
% $$
%
% is referred to as "-Psi d2F/dpdp-" and currently implemented for $p=C$.
%
% The term
%
% $$
%   \frac{\partial^2 J}{\partial q_k\, \partial q_m} \xi_{ki} \, \xi_{mj}
% $$
%
% is referred to as "-xi d2J/dqdq xi-" and implemented for $A$, $B$, and $C$.
%
%
% The term
%
% $$
%   \Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial q_m} \xi_{ki} \, \xi_{mj}
% $$
%
% is not implemented.
%
% The term
%
% $$
%   \frac{\partial^2 J}{\partial p_i \, \partial q_k} \xi_{kj}
% $$
%
% is zero for
%  $J=R(p) + I(q(p))$,
% (these are partial derivatives).
%
% The term
%
% $$
%   \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj}
% $$
%
% is not implemented (but looks relatively easy).
%
% The term
%
% $$
%   \frac{\partial^2 J}{\partial q_k \, \partial p_j} \xi_{ki}
% $$
%
% is zero for
%  $J=R(p) + I(q(p))$,
% (these are partial derivatives).
%
% The term
%
% $$
%   \Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki}
% $$
%
% is not implemented (but similar to an above term, and looks relatively easy.)
%
%
%  Currently only for u, v as q variables
%
%
% $$2 I= u M_u u + v M_v v + \dot{h} M_{\dot{h}} \dot{h} $$
%
% $$ p =\left ( \begin{array}{c}  A \\ B \\ C  \end{array} \right ) $$
%
% $$ q =\left ( \begin{array}{c}  u \\ v \\ \dot{h}  \end{array} \right ) $$
%
% $$ \xi = \frac{\partial q}{\partial p} = \left ( \begin{array}{cc}
%  \frac{\partial u}{ \partial A} & \frac{\partial u}{\partial C}  \\
%  \frac{\partial v}{\partial A } & \frac{\partial v}{\partial C}  \\
%  \frac{\partial \dot{h}}{\partial A} & \frac{\partial \dot{h}}{\partial C}  \\
% \end{array} \right ) $$
%
% $$ \frac{\partial^2 J}{\partial q\, \partial q} = \left ( \begin{array}{ccc}
%                                                    M_u & 0 & 0 \\
%                                                    0 & M_v & 0 \\
%                                                    0 & 0 & M_{\dot{h}}
% \end{array} \right ) $$
%
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
%
%
% https://www.sciencedirect.com/science/article/pii/S0377042708006523
%
%%

%% Get the sensitivity matrices



%CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;

H=0 ;
HessianTerms="-xi d2J/dqdq xi-" ;   % this often results in amazingly good convergence! Many order of magnitude decrease per iteration, for example in the AC inversion test 
                                    % the cost function goes from 1e5 to 1e-25 in 4 iterations
% HessianTerms="-Psi d2F/dpdp-";    % While this does work, the performance is not particularly good, maybe 50% reduction per
                                    % iteration, and often not much better than the gradient descent. 
% HessianTerms="-xi d2J/dqdq xi-Psi d2F/dpdp-"; % when adding "-Psi d2Fdpdp-" the inversion performs worse...! I suspect
                                                 % that something is not quite right with the -Pis d2F/dpdp-" calculations 


if contains(HessianTerms,"-xi d2J/dqdq xi-")
    GetSensitivites=true;
else
    GetSensitivites=false;
end




if GetSensitivites

    if  ~contains(CtrlVar.Inverse.Measurements,'-uv-','IgnoreCase',true)

        fprintf("Currently, when calculating the Hessian using the Direct-Adjoint method, measurement must include surface velocities.\ n")
        fprintf("One can either only use velocity measurements (uv), or measurements of both velocity and rate of thickness changes (uv and dhdt), but not thickness changes alone. \n")
        error("CalcDirectAdjointHessian:noVelMeas","Measurements must include surface velocities")


    end


  
    %
    % if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    %     [KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC]=duvhdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs) ;
    % else

    [KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC]=duv_hdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs);
    % end

end


if contains(HessianTerms,"-xi d2J/dqdq xi-")

    % This requires too large memory, possible approach is to calculate Hessian-vector product and only one row of the Hessian at
    % a time.
    %
    % https://arxiv.org/html/2410.22575v1
    %
    %
    %

    if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)

        xi=[KdudA KdudB KdudC ; KdvdA KdvdB KdvdC ; KdhdA KdhdB KdhdC] ; % 2 Nnodes \times nP Nnodes where nP is the number of fields inverted for, e.g. 2 if inverting for A and C, 1 if only inverting for B

        d2Jdqq=blkdiag(d2Iduu,d2Idvv,d2Idhdothdot);  % d2Iduv and d2Idvu are zeros, 2 Nnodes \times 2 Nnodes

    else

        xi=[KdudA KdudB KdudC ; KdvdA KdvdB KdvdC] ; % 2 Nnodes \times nP Nnodes where nP is the number of fields inverted for, e.g. 2 if inverting for A and C, 1 if only inverting for B

        d2Jdqq=blkdiag(d2Iduu,d2Idvv);  % d2Iduv and d2Idvu are zeros, 2 Nnodes \times 2 Nnodes

    end

%%
%
% $$ 
% \xi^T \frac{\partial^2 J }{\partial q \, \partial q} \, \xi
% $$
%
%
%
%%

    tMult=tic;
    H=xi'*(d2Jdqq*xi)+H ;

    H=0.5*(H+H');
    tMult=toc(tMult);
    fprintf(" Multiplication calculated in %f sec\n",tMult)

end

if contains(HessianTerms,"-Psi d2F/dpdp-")

    H=H+PsiTimesddFuvdpdp(CtrlVar,MUA,F,uAdjoint,vAdjoint);

end


end