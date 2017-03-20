function [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,plb,pub,RunInfo)


CtrlVar.Inverse.MatlabOptimisationParameters= optimoptions(CtrlVar.Inverse.MatlabOptimisationParameters,'MaxIterations',CtrlVar.Inverse.Iterations);


Test=CtrlVar.Inverse.MatlabOptimisationParameters;

if isa(Test,'optim.options.Fminunc')
    
    [p,J,exitflag,output] = fminunc(func,p0,CtrlVar.Inverse.MatlabOptimisationParameters);
    
    if isfield(RunInfo.Inverse,'fminunc')
        RunInfo.Inverse.fminunc=output;
    end
    
elseif isa(Test,'optim.options.Fmincon')
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    
    
    [p,J,exitflag,output] = fmincon(func,p0,A,b,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationParameters);
    
    if isfield(RunInfo.Inverse,'fmincon')
        RunInfo.Inverse.fmincon=output;
    end
    
else
    
    fprintf('Matlab Optimization selected, but Matlab optimization routine not recognized.\n')
    fprintf(' Either select fminunc or fmincon. \n')
    error(' invalid input parameters ')
    
end

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