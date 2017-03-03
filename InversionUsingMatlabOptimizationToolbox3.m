function [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,RunInfo)



[p,J,exitflag,output] = fminunc(func,p0,CtrlVar.Inverse.MatlabOptimisationParameters);

[stop,Outs] = fminuncOutfun();


RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+Outs.iteration];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Outs.fval];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.J;Outs.StepSize];
RunInfo.Inverse.R=[RunInfo.Inverse.R;Outs.fval+NaN];
RunInfo.Inverse.I=[RunInfo.Inverse.I;Outs.fval+NaN];
RunInfo.Inverse.p=Outs.p;
% If I need some further info and want to update F
%[J,Gradient,Hessian,Outs,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);



end