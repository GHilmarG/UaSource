





function H = CalcDirectAdjointHessian(UserVar,CtrlVar,RunInfo,MUA,F,BCs,l,Meas,BCsAdjoint,Psi_x,Psi_y)

narginchk(11,11)


%% Calculates the Hessian, H, using the direct-adjoint approach.
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
% is easy, and is done in Regularisation.m
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
% $$ F^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi_i \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle t_x | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% $$ F^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi_i \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle t_y | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
% Here we use
%
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
% For Weertman:
%
% $$t_x=\mathcal{G} \beta^2  u$$
%
% $$t_y=\mathcal{G} \beta^2  v$$
%
% $$ \beta^2 = (C+C_0)^{-1/m} \; (u^2+v^2+v_0^2)^{(1/m-1)/2} $$
%
%
% *Terms:*
%
% The whole expression is:
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
% *First term:*
%
% $$ \frac{\partial^2 J}{\partial p_i \, \partial p_j} $$
%
% This is simple and only involves explicit dependency on p. In the case of $A$, $B$ and $C$ this involves the precision matrix of the prior.
% And in the case of $B$ there is an additional contribution from the direct measurements misfit term and this is done in the
% misfit function.
%
% *Second term:*
%
% $$ \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial p_j} $$
%
% This is done in separate m-files for A, B and C. It is zero if the problem is linear in $p$. The forward problem is not
% linear in neither $A$ or $C$, even for $n=1$ and $m=1$.
%
% *Third term:*
%
% $$ \frac{\partial^2 J}{\partial q_k\, \partial q_m} \xi_{ki} \, \xi_{mj} $$
%
% This involves second-order derivatives of $J$ with respect to A and C which are the misfit precision matrices. This is
% (almost) always a very important term and should be included. It is currently implemented for $A$, $B$ and $C$.
%
% *Forth term*
%
% $$ \Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial q_m} \xi_{ki} \, \xi_{mj} $$
%
% This involves second-order (explicit) derivatives of the forward model with respect to u,v, and h.
%
% *Fifth term*
%
% $$ \frac{\partial^2 J}{\partial p_i \, \partial q_k} \xi_{kj} $$
%
% This involves a mixed derivative of the cost function. The $J$ terms typically either involve only $p$ or $q$ but not both.
% For example $C$ will be involved in the prior and $u$ in the misfit term, but those are separate terms. Therefore this term
% in the Hessian is zero.
%
% *Sixth term*
%
% $$  \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj} $$
%
% Here we have mixed derivatives of the forward model. These are not zero, even if the problem is linear in both $p$ and $q$.
%
%
% *Seventh term*
%
% $$\frac{\partial^2 J}{\partial q_k \, \partial p_j} \xi_{ki} $$
%
% Mixed derivatives of the const function. As above, this term is zero.
%
% *eighth term*
%
% $$ \Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki} $$
%
% This term also involves mixed derivatives of the forward model, and is not zero
%
% If
%
% $$ H^{\#8}_{ij}=\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki} $$
%
% and
%
% $$ H^{\#6}_{ij} = \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj} $$
%
% then
%
% $$ [H^{\#8}]^T=H^{\#8}_{ji}=\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_i} \xi_{kj} =\Psi_n \frac{\partial^2 F_n}{\partial p_i \,\partial q_k} \xi_{kj} = H^{\#6}_{ij}= [H^{\#6}] $$
%
%
%
%
%%

%% Get the sensitivity matrices



%CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;

H=0 ;

HessianTerms="-xi Jqq xi-xi Fqq xi-Fpp-Fpq xi-" ;

if contains(HessianTerms,"-xi Jqq xi-") || contains(HessianTerms,"-xi Fqq xi-")
    GetSensitivites=true;
else
    GetSensitivites=false;
end



%% sensitivity matrix, \xi = \partial q / \partial p   % tested
if GetSensitivites

    if  ~contains(CtrlVar.Inverse.Measurements,'-uv-','IgnoreCase',true)

        fprintf("Currently, when calculating the Hessian using the Direct-Adjoint method, measurement must include surface velocities.\ n")
        fprintf("One can either only use velocity measurements (uv), or measurements of both velocity and rate of thickness changes (uv and dhdt), but not thickness changes alone. \n")
        error("CalcDirectAdjointHessian:noVelMeas","Measurements must include surface velocities")


    end

    [KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC]=duv_hdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs);
    xi=[KdudA KdudB KdudC ; KdvdA KdvdB KdvdC] ;



    %% H^{qq}
    %
    % $$\xi^T (J^{qq}+\mathcal{F}^{qq} )\xi$$
    %

    KJqq=0; KFqq=0; 

    if contains(HessianTerms,"-xi Jqq xi-")

        KJqq=Jqq(CtrlVar,MUA,F,BCs,Meas);

    end


    if contains(HessianTerms,"-xi Fqq xi-")

        KFqq=Fqq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
    end

    KJqqFqq=KJqq+KFqq;


    % the numerical sparsity of xi is close to 1 (ie not sparse at all)
    % much better to use full matrix in the multiplication
    %
    % requires
    %

    % tMult=tic;
    % H=xi'*(d2Jdqq*xi)+H ;
    %
    % H=0.5*(H+H');
    % tMult=toc(tMult);
    % fprintf(" Multiplication calculated in %f sec\n",tMult)

    xiNumericalSparsity=nnz(xi)/numel(xi);
    if xiNumericalSparsity>0.5
        xi=full(xi);
    end


    tMult=tic;
    H=xi'*(KJqqFqq*xi)+H ;


    tMult=toc(tMult);
    %  fprintf(" Multiplication calculated in %f sec\n",tMult)
end

%%  H^{pp}  (but only the F^pp contribution)  : Tested
%
% $$H^{p}=J^{pp}+\mathcal{F}^{pp}$$
%
%
if contains(HessianTerms,"-Fpp-")  % this is from $\delta^2_{pp} F$

    KFpp=Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
    H=H+KFpp; % Here missing is the Jpp contribution, but this is added in the Regularisation.m function

end


%% H^{pq}+H^{qp}  : Tested
%
% $$H^{pq}+H^{qp}=\big(J^{pq}+\mathcal{F}^{pq}\big)\xi \;+\; \Big[\big(J^{pq}+\mathcal{F}^{pq}\big)\xi\Big]^T$$
%
% The J^{pq} contribution is missing, but this is likely to be zero as each term in the cost function is only an explicit
% function of either p or q, not both.
%
if contains(HessianTerms,"-Fpq xi-") % this is from $\delta^2_{pq} F$ and $\delta^2_{qp} F $


    [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC);
   
    H=H+KHess_qp ;

end


H=0.5*(H+H');

end