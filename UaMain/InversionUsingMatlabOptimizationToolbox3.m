




function   [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc,Aineq,bineq)



CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);


CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);

if CtrlVar.Inverse.MatlabOptimisationHessianParameters.Algorithm=="trust-region-reflective"

    % I don't think this is needed when using the trust-region-reflective, because this algorithm uses the third
    % argument of @func which here is JGH

    if CtrlVar.Inverse.MinimisationMethod=="-MatlabOptimization-DirectAdjointHessian-"


        %% Hessian provided

        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn','objective');
        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianMultiplyFcn',[]);

    elseif CtrlVar.Inverse.MinimisationMethod=="-MatlabOptimization-HessianVectorProduct-"

        %% Hessian-vector product

        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',[]);
        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianMultiplyFcn', @(Hinfo,v) HessianVectorProduct(Hinfo,v,func));

    elseif CtrlVar.Inverse.MinimisationMethod=="-MatlabOptimization-HessianFiniteDifferences-"

        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',[]);
        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianMultiplyFcn',[]);
        
       %HessPattern=speye(numel(p0),numel(p0));    % diagonal
       % HessPattern=spdiags([1 1 1],-1:1,numel(p0),numel(p0));  % tri-diagonal 

       % As expected, the convergence does improve as the Hessian sparsity is decreased (higher n), however, not significantly so
       % and the convergence was always linear. 
        n=11; 
        HessPattern=spdiags(ones(1,n),-(n-1)/2:(n-1)/2,numel(p0),numel(p0));  % n-diagonal 
        
        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessPattern',HessPattern);
        CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FiniteDifferenceType','central');
        % v=zeros(numel(p0),1)+1e100;
        % CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FiniteDifferenceStepSize',v);
        


    end

else
    CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',Hfunc);
end


Test=CtrlVar.Inverse.MatlabOptimisationGradientParameters;




if isa(Test,'optim.options.Fminunc')

    [p,J,exitflag,output] = fminunc(func,p0,CtrlVar.Inverse.MatlabOptimisationGradientParameters);

    if isfield(RunInfo.Inverse,'fminunc')
        RunInfo.Inverse.fminunc=output;
    end

elseif isa(Test,'optim.options.Fmincon')


    Aeq = [];
    beq = [];
    nonlcon = [];


    if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")



        [p,J,exitflag,output,lambda,grad,hessian] = fmincon(func,p0,Aineq,bineq,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationHessianParameters);

    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Gradient")


        [p,J,exitflag,output] = fmincon(func,p0,Aineq,bineq,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationGradientParameters);

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