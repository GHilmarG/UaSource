




function [J,dJdp,Hessian,JGHouts,F]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)



%%
%
% JGH: Returns the cost function (J), the gradient of the cost function with respect to p (dJdp), and the Hessian (ddJddp).
%
% The Hessian of the regularization term (R) can usually be calculated exactly, while the Hessian of the misfit/likelihood
% term (I), can not. However, one can come up with a educated guess for the Hessian of I with respect to C.
%
%
% Calculates objective function (J), gradient (dJdp, accurate), Hessian (guessed).
%
%
%%

persistent ubP vbP JGH1 JGH2 JGH3

narginchk(11,11)


%% some counters for how often JGH is called and with what number of arguments
if isempty(JGH1)
    JGH1=0;  % counter for 1-argument output, just cost function evaluation
    JGH2=0;  % counter for 2-argument output, cost and gradient
    JGH3=0;  % counter for 3-argument output, cost, gradient and Hessian
end


JGH1=JGH1+1;


%% Create some internal flags indicating if cost, gradient or Hessian need to be calculated and returned.

CtrlVar.Inverse.CalcGrad=false;
CtrlVar.Inverse.CalcGradI=false;
CtrlVar.Inverse.CalcGradR=false;


CtrlVar.Inverse.CalcHess=false;
CtrlVar.Inverse.CalcHessI=false;
CtrlVar.Inverse.CalcHessR=false;

if nargout==1
    dJdp=[] ; Hessian=[] ; JGHouts=[] ;
end

if nargout>=2  % always calculates the gradient if number of output arguments is 2 or larger
    CtrlVar.Inverse.CalcGrad=true;
    CtrlVar.Inverse.CalcGradI=true;
    CtrlVar.Inverse.CalcGradR=true;
end

% Calculate the Hessian provided:
%     1) CtrlVar.JGH.CalcHessian=true;
% and 2) number of output arguments is larger or equal to 3
if nargout>= 3  &&  CtrlVar.JGH.CalcHessian
    CtrlVar.Inverse.CalcHess=true;
    CtrlVar.Inverse.CalcHessI=false;
    CtrlVar.Inverse.CalcHessR=false;
else
    CtrlVar.Inverse.CalcHess=false;
    CtrlVar.Inverse.CalcHessI=false;
    CtrlVar.Inverse.CalcHessR=false;
    Hessian=[];
end


%%

[isA,isB,isC] = isABC(CtrlVar);
[is_uv_meas,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar);



%% The function requires a solution of the forward model, often with just a slightly different input parameters.
% Therefore, save previous velocity solution as persistent variables and use as initial staring point for next uv-solve.

if CtrlVar.Inverse.ResetPersistentVariables
    ubP=[];
    vbP=[];
end

if ~isempty(ubP)
    F.ub=ubP;
    F.vb=vbP;
end




% The vector p contains the variables for which the inversion is being performed. So if the inversion is done over log(c)
% only, then p=log(C). And if the inversion is done over A, B and C then p=[A;B;C].

% Populate F with the current values in the vector p ahead of a call the the forward model.
F=p2F(CtrlVar,MUA,p,F,Meas,Priors);  % This maps from the vector p to the field variables F
% p is the vector of control variables. This is some combination of A, B and C
% When inverting for A, we have p=[log10(AGlen]
% When inverting for C, we have p=[log10(C)]
% When inverting for A and C we have p=[log10(AGlen);log10(C)]
% and so on


%% Forward model solution
[~,~,F,l,dFduv]= uv([],[],CtrlVar,MUA,BCs,F,l);

if is_dhdt_meas
    % If dh/dt is included as measurements, I need that calculated dh/dt, which in turn requires the mass-balance, a, as well
    if isempty(F.as) || isempty(F.ab)
        [~,F]=GetMassBalance([],CtrlVar,MUA,F);
    end
    [~,F.dhdt]=dhdtExplicit([],CtrlVar,MUA,F,BCs) ;
end

%% Cost function and gradient
% The cost function, J), is split into a misfit (I) and a regularization term (R). These usually consist of further
% terms.
%
% Get the I and R terms, and the gradients if required. The Hessian is calculated later.


[R,dRdp]=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);
[I,dIdp,Psi_x,Psi_y,F,dFduv]=Misfit(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,dFduv);

if  CtrlVar.Inverse.CalcGrad  % gradient needed
    dJdp=dRdp+dIdp;
    JGH2=JGH2+1;

    %% Box transformation: chain rule
    %
    % dJ/du = (dp/du) * dJ/dp , with dp/du = (ub-lb)*r*(1-r) , see BoxTransform.m .
    %
    % This is a relation between l2 gradients, so it MUST be applied here, to dRdp+dIdp, and before the Riesz map below.
    % Applied after the Riesz map it would read G\(D.*(G*dJdp)) instead of D.*dJdp, costing an extra matrix-vector
    % product and an extra solve with G on every cost-function evaluation.
    %
    % The Jacobian is evaluated at u. The vector p arriving here is the variable the optimiser holds, which is u if only
    % the box transformation is used, and the Cholesky-mapped u if that mapping is used as well. In the latter case u has
    % to be recovered first, by the same inverse mapping p2F applies.

    if CtrlVar.Inverse.BoxTransform

        if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient
            u=MUA.PRG*(MUA.RG\p) ;
        else
            u=p ;
        end

        dJdp=BoxTransform("jacobian",u,MUA.BoxTransform).*dJdp ;

    end

    if CtrlVar.Inverse.RieszMapGradient

        if ~isfield(MUA,"dG") && isempty(MUA.dG)
            dJdp=MUA.G\dJdp; 
        else
            dJdp=MUA.dG\dJdp;
        end

        % To do: Here I could add the BCs in the future

    end

    if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient
       
        %fprintf("Cholesky-mapped gradient.\n")
     
        %dJdp=MUA.RG*dJdp; 
        dJdp=MUA.RG*(MUA.PRG'*dJdp);

    end

end

%% Hessian
if CtrlVar.Inverse.CalcHess  % Hessian needed

    Hessian = BuildInversionHessian(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y);
    assert(numel(dJdp)==size(Hessian,1),"Regularisation:DimentionalMismatch","sizes of gradient and Hessian not compatible.")
    JGH3=JGH3+1;
end


if F.solution=="-uv-"
    % To speed up the forward solve, the previous solution is stored locally and then used as a starting value in next
    % calculation. The idea is that usually the parameter vector (p) only changes slightly form one inverse iteration to the
    % next, so the (u,v) solution is likely to be similar to the previously calculated one.
    ubP=F.ub;
    vbP=F.vb;
else
    warning('JGH:returnsNaN',' uv solution did not converge. Returning NaN in cost function.\n ') ;
    ubP=[];
    vbP=[];
    I=NaN;
    R=NaN ;
    dJdp=p*0+NaN;
    
end

J=full(R+I);

if J < 0
    fprintf("J less than zero!! \n")
end


if nargout>3  % additional information needed as output
   
    JGHouts.J=J;
    JGHouts.R=R;
    JGHouts.I=I;
 
    JGHouts.dRdp=dRdp;
    JGHouts.dIdp=dIdp;
    JGHouts.ddIdpp=[];
    JGHouts.ddRdpp=[];
    JGHouts.Psi_x=Psi_x;
    JGHouts.Psi_y=Psi_y;

else
    JGHouts=[];
end



if isnan(J)
    warning("JGH:ObjectivFunctionIsNaN","objective function is nan")
end


% fprintf("JGH(%i,%i,%i): \t J=%g\n",JGH1,JGH2,JGH3,J)

end

