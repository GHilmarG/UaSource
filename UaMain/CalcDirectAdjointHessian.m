





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
% Limitations: Currently the Hessian calculations have been implemented for logA and logC, but not for B
%
%
%%


% [~,~,F,l]= uv([],[],CtrlVar,MUA,BCs,F,l);

%%


if CtrlVar.Inverse.BoxTransform
    fprintf("CtrlVar.Inverse.BoxTransform=true, but this is not yet implemented for the Direct-Adjoint Hessian approach.\n")
    fprintf("CalcDirectAdjointHessian:Not implemented for box transform.\n")
    error("NotImplemented")
end
%% I label individual Hessian terms and have the option of only calculating a subset of those for testing purposes.

% HessianTerms="-xi Jqq xi-xi Fqq xi-Fpp-Fpq xi-Jpp-" ;

HessianTerms=CtrlVar.Inverse.HessianTerms;

if contains(lower(HessianTerms),"-all-")
    HessianTerms="-xi Jqq xi-xi Fqq xi-Fpp-Fpq xi-Jpp-" ;
end

%% Do I need to calculate the sensitivity matrices?


H=[] ;

% Do I need to calculate the sensitivities, i.e. the Jacobian dq/dp, for the Hessian terms requested? 
if contains(HessianTerms,"-xi Jqq xi-") || contains(HessianTerms,"-xi Fqq xi-")
    Sensitivites=true;
else
    Sensitivites=false;
end



%% sensitivity matrix, \xi = \partial q / \partial p   % tested
if Sensitivites
    % Note: these have been calculated for logA, B, and logC. So this can be considered done. 
    [KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC]=duv_hdABC(CtrlVar,MUA,F,l,BCs);
    xi=[KdudA KdudB KdudC ; KdvdA KdvdB KdvdC] ;
else
    xi=[];

end

%% H^{qq}
%
% $$\xi^T (J^{qq}+\mathcal{F}^{qq} )\xi$$
%

KJqq=0; KFqq=0;

if Sensitivites

    if contains(HessianTerms,"-xi Jqq xi-")
        % Note: These are done for all measurement types, u,v and dot{h}. So this is done, although the expressions are not quite
        % correct for spatially variable measurement errors. 
        %
        % These are explicit second-order derivatives of the cost function, J, with respect to the q=(u,v)
        %
        % There is a term related to dot(h) but this is really due to that term being an explicit function of q=(u,v)
        %
        KJqq=Jqq(CtrlVar,MUA,F,BCs,Meas);
    end

    if contains(HessianTerms,"-xi Fqq xi-")
        %
        % Note: These are explicit second-order derivatives of the forward model with respect to q=(u,v)
        %
        % 
        %
        KFqq=Fqq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
    end
    KJqqFqq=KJqq+KFqq;

    xiNumericalSparsity=nnz(xi)/numel(xi);
    if xiNumericalSparsity>0.5
        xi=full(xi);
    end


    H=xi'*(KJqqFqq*xi) ;

end
%  fprintf(" Multiplication calculated in %f sec\n",tMult)


%%  H^{pp} , )
%
% $$H^{p}=J^{pp}+\mathcal{F}^{pp}$$
%
% F^pp contribution
if contains(HessianTerms,"-Fpp-")  % this is from $\delta^2_{pp} F$


    % Note: These are second-order derivatives of the forward model with respect to p=(logA,B,logC). 
    % This has not been implemented for B
    KFpp=Fpp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;
    if isempty(H)
        H=KFpp;
    else
        H=H+KFpp;
    end
end

% Jpp
if contains(HessianTerms,"-Jpp-")  % explicit dependency of J on p=(logA,B,logC)

    % Note: these are second-order derivatives of the cost function with respect to p=(logA,B,logC)
    % This has been implemented for all three p fields
    KJpp=Jpp(CtrlVar,MUA);

    if isempty(H)
        H=KJpp;
    else
        H=H+KJpp;
    end


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

if Sensitivites
    if contains(HessianTerms,"-Fpq xi-") % this is from $\delta^2_{pq} F$ and $\delta^2_{qp} F $

        % These are second-order mixed derivatives of the forward model with respect to q,p
        %
        % This has not been implemented for derivatives involving B
        [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC);
        H=H+KHess_qp ;

    end
end

H=0.5*(H+H');



if CtrlVar.Inverse.TestDirectAdjoint.isTrue
    FiniteDifferenceTestAndPlots(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,H)
end


end


function   FiniteDifferenceTestAndPlots(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,H)


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


figDA=FindOrCreateFigure("Test: Direct-Adjoint H") ; clf(figDA)


plot(Hcolumn,H_FD,"or") ; axis equal ;
hold on ;
plot([min(Hcolumn) max(Hcolumn)],[min(Hcolumn) max(Hcolumn)],"--k")

ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$H$, Direct-Adjoint approach",Interpreter="latex")
subtitle(sprintf("Comparison is here for one random column: %i",iColumn),Interpreter="latex")



end


