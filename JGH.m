function [J,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)




% Calculates objective function, gradient (accurate), Hessian (guessed)

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


if contains(CtrlVar.Inverse.MinimisationMethod,"UaOptimization")
    % p=kk_proj(p,pub,plb);  % I guess the matlab optimisation toolbox uses a bit more sophisticated approach (I hope).
    %[~,iU,iL] = kk_proj(p,pub,plb);
    %numel(find(iU))
    %numel(find(iL))
end

% Note: I should consider writing this as F=p2InvValues(CtrlVar,p)

F=p2F(CtrlVar,MUA,p,F,Meas,Priors);

if any(isnan(F.C)) 
    save TestSave ; 
    error( ' C nan ') ; 
end

[UserVar,RunInfo,F,l,dFduv]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);


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

