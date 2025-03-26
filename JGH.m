




function [J,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)


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

persistent ubP vbP

narginchk(14,14)
CtrlVar.nargoutJGH=nargout;

if nargout==1
    CtrlVar.Inverse.CalcGradI=false;
    CtrlVar.Inverse.CalcGradR=false;
    dJdp=[] ; Hessian=[] ; JGHouts=[] ;
else
    CtrlVar.Inverse.CalcGradI=true;
    CtrlVar.Inverse.CalcGradR=true;
end


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
F=p2F(CtrlVar,MUA,p,F,Meas,Priors);

if any(isnan(F.C)) 
    save TestSave ; 
    error( ' C nan ') ; 
end

% This is a call to the forward model.
[UserVar,RunInfo,F,l,dFduv]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);


% The cost function, J), is split into a misfit (I) and a regularization term (R). These usually consist of further
% terms.
%
% Get the I and R terms, and the gradients if required.
if nargout==1
    R=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo) ;
    I=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dFduv) ;
else
    [R,dRdp,ddRddp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo) ;
    [I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dFduv) ;
end



if nargout>1   % gradient needed
    dJdp=dRdp+dIdp;
end

if nargout>2  % Hessian needed 
    if isempty(ddIddp)
        Hessian=ddRddp;
    else
        Hessian=ddRddp+ddIddp;
    end
end


if RunInfo.Forward.uvConverged
     % To speed up the forward solve, the previous solution is stored locally and then used as a starting value in next
     % calculation. The idea is that usually the parameter vector (p) only changes slightly form one inverse iteration to the
     % next, so the (u,v) solution is likely to be similar to the previously calculated one.
    ubP=F.ub; 
    vbP=F.vb;
else
    warning('JGH:returninNaN',' uv solution did not converge. Returning NaN in cost function.\n ') ;
    ubP=[];
    vbP=[];
    I=NaN;
    R=NaN ;
    dJdp=p*0+NaN;
    MisfitOuts.I=NaN;
end

J=R+I;

if nargout>3
    JGHouts.dRdp=dRdp;
    JGHouts.dIdp=dIdp;
    JGHouts.ddIdpp=ddIddp;
    JGHouts.ddRdpp=ddRddp;
    JGHouts.RegOuts=RegOuts;
    JGHouts.MisfitOuts=MisfitOuts;
else
    JGHouts=[];
end



end

