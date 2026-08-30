





function H = CalcDirectAdjointHessian(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y)



narginchk(10,10)


%% Calculates the Hessian, H, using the direct-adjoint approach.
%
%
% $$
% H_{lm}  = \delta^2_{qq} J[\xi_{,l},\xi_{,m}] 
% + \delta^2_{qp}J[\xi_{,l},\phi_m] 
% + \delta^2_{pq}J[\phi_l,\xi_{,m}] 
% + \delta^2_{pp}J[\phi_l,\phi_m]  
%  + \langle \Psi ,
% \delta^2_{qq}\mathcal{F}[\xi_{,l},\xi_{,m}] 
%     + \delta^2_{qp}\mathcal{F}[\xi_{,l},\phi_m] 
%     + \delta^2_{pq}\mathcal{F}[\phi_l,\xi_{,m}] 
%     + \delta^2_{pp}\mathcal{F}[\phi_l,\phi_m] 
%     \rangle
% $$
%
% These terms can be grouped together and renames as:
%
%
%
% $$
% H = \underbrace{\xi^T\big(J^{qq}+\mathcal{F}^{qq}\big)\xi}_{H^{qq}}
% \;+\; \underbrace{\big(J^{pq}+\mathcal{F}^{pq}\big)\xi \;+\; \Big[\big(J^{pq}+\mathcal{F}^{pq}\big)\xi\Big]^T}_{H^{pq}+H^{qp}}
% \;+\; \underbrace{\big(J^{pp}+\mathcal{F}^{pp}\big)}_{H^{pp}}
% $$
%
%
% If the forward model is
%
% $$ \mathcal{F}(q(p),p) = 0 $$
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
%
%%


% [~,~,F,l]= uv([],[],CtrlVar,MUA,BCs,F,l);



%% I label individual Hessian terms and have the option of only calculating a subset of those for testing purposes.

HessianTerms="-xi Jqq xi-xi Fqq xi-Fpp-Fpq xi-" ;

%% Do I need to calculate the sensitivity matrices?


H=0 ;


if contains(HessianTerms,"-xi Jqq xi-") || contains(HessianTerms,"-xi Fqq xi-")
    GetSensitivites=true;
else
    GetSensitivites=false;
end



%% sensitivity matrix, \xi = \partial q / \partial p   % tested
if GetSensitivites


    [KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC]=duv_hdABC(CtrlVar,MUA,F,l,BCs);
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

%%  H^{pp} , (here only the F^pp contribution)  : Tested
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
% The J^{pq} contribution is not missing as each term in the cost function is only an explicit
% function of either p or q, not both. 
% 
% Even the $$J_{\dot{h}}$$ terms only involves $u$ and $v$ and not any of $A$, $B$ or
% $C$, so here $$J_{\dot{h}}^{pq} =0 $$ as well
%
if contains(HessianTerms,"-Fpq xi-") % this is from $\delta^2_{pq} F$ and $\delta^2_{qp} F $


    [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC);
   
    H=H+KHess_qp ;

end


H=0.5*(H+H');


% %% J^{pp} : I still have the J^pp to add, but this is done in the Regularisation part and added later in the code. But I
% should consider changing this and make sure all contributions are included here. The missing contribution is ddRdpp as
% returned by: 
% 
%  [R,dRdp,ddRdpp]=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ; 
% 
% 
% %%



if CtrlVar.Inverse.TestDirectAdjoint.isTrue

    FiniteDifferenceTestAndPlots(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,H)

end



end


function   FiniteDifferenceTestAndPlots(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,H)


 % Since the above calculation of the Hessian does not include Jpp, I better add that here
[R,dRdp,ddRdpp]=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ; 

H=H+ddRdpp;

% First map all A and C fields to p. This takes care of the log conversion
[p,plb,pub]=F2p(CtrlVar,MUA,F); 

% the do the perturbation with respect to p

iColumn=randi(numel(p));
%iColumn=1209; 

% Perform perturbation on the selected column
perturbation = 1e-3; % Define a small perturbation value. Be careful that this is in log space for A and C. Might need to try out several different amplitudes

pPerturbed_pos = p; 
pPerturbed_pos(iColumn) = pPerturbed_pos(iColumn) + perturbation;


% I now map to F from p. However, this should not be needed as this is always done in JGH. The reason JGH needs F as an input
% at all is because F contains various other fields that are not dependent on p, but I still need those as input fields for
% the forward model.

% F=p2F(CtrlVar,MUA,pPerturbed_pos,F,Meas,Priors); 

% JGH calculates the cost function (J), the gradient (G) and the Hessian (H). Here I only need the gradient.
% 
% Note: If I were to include a third output argument, which is the Hessian, the JGH function would call
% CalcDirectAdjointHessian.m, resulting in an endless recursion.

[J_pos,dJdp_pos]=JGH(pPerturbed_pos,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);

pPerturbed_neg = p; 
pPerturbed_neg(iColumn) = pPerturbed_neg(iColumn) - perturbation;

% F=p2F(CtrlVar,MUA,pPerturbed_neg,F,Meas,Priors); 

[J_neg,dJdp_neg]=JGH(pPerturbed_neg,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);

H_FD=(dJdp_pos-dJdp_neg)/(2*perturbation) ;

Hcolumn=H(:,iColumn);

Diff=norm(Hcolumn-H_FD)/norm(Hcolumn);
fprintf("H: normalized norm of difference between Direct-Adjoint and FD for column %i is %g \n",iColumn,Diff)


figDA=FindOrCreateFigure("Test: Direct-Adjoint") ; clf(figDA)


plot(Hcolumn,H_FD,"or") ; axis equal ;
hold on ;
plot([min(Hcolumn) max(Hcolumn)],[min(Hcolumn) max(Hcolumn)],"--k")

ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$H$",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")



end


