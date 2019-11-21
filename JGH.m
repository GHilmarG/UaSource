function [J,dJdp,Hessian,JGHouts,F]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

% Calculates objective function, gradient (accurate), Hessian (guessed)

persistent ubP vbP

narginchk(14,14)
CtrlVar.nargoutJGH=nargout;

if nargout==1
    CtrlVar.Inverse.CalcGradI=false;
    CtrlVar.Inverse.CalcGradR=false;
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


if CtrlVar.Inverse.MinimisationMethod=="UaOptimization"
     p=kk_proj(p,pub,plb);  % I guess the matlab optimisation toolbox uses a bit more sophisticaed approach (I hope). 
end

% Note: I should consider writing this as F=p2InvValues(CtrlVar,p)

F=p2F(CtrlVar,MUA,p,F,Meas,Priors); 



[UserVar,RunInfo,F,l,dFduv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

[R,dRdp,ddRddp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo) ;
[I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dFduv) ;



if nargout>1
    dJdp=dRdp+dIdp;
    Hessian=ddRddp+ddIddp;
end

if RunInfo.Forward.Converged
    ubP=F.ub;
    vbP=F.vb;
else
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

