
function   [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc)

CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);


CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);



CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',Hfunc);

Test=CtrlVar.Inverse.MatlabOptimisationGradientParameters;

CtrlVar.Inverse.MinimisationMethod="MatlabOptimization-GradientBased";  % rhubarb 


if isa(Test,'optim.options.Fminunc')
    
    [p,J,exitflag,output] = fminunc(func,p0,CtrlVar.Inverse.MatlabOptimisationGradientParameters);
    
    if isfield(RunInfo.Inverse,'fminunc')
        RunInfo.Inverse.fminunc=output;
    end
    
elseif isa(Test,'optim.options.Fmincon')
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    
    
    if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
  
        
        
        [p,J,exitflag,output] = fmincon(func,p0,A,b,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationHessianParameters);
        
    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Gradient")
        
        
        [p,J,exitflag,output] = fmincon(func,p0,A,b,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationGradientParameters);
        
    else

        fprintf("The variable CtrlVar.Inverse.MinimisationMethod has an invalid value. ") 
        error("InversionUsingMatlabOptimizationToolbox3:InvalidParameters","CtrlVar.Inverse.MinimisationMethod invalid.")
        
    end
    
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
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;Outs.GradNorm];
RunInfo.Inverse.p=Outs.p;
% If I need some further info and want to update F
%[J,Gradient,Hessian,Outs,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);



end