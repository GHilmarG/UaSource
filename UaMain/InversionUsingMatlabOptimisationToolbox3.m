




function   [p,RunInfo]=InversionUsingMatlabOptimisationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc,Aineq,bineq)


if CtrlVar.Inverse.RieszMapGradient 
    error("InversionUsingMatlabOptimisationToolbox3:NoRieszMappingAllowedWithMATLABfmincon","Do not compbine Riesz mapping and fmincon. But you can combine Cholesky mapping with fmincon.")
end




CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationGradientParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationGradientParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);


CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'MaxIterations',CtrlVar.Inverse.Iterations);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'OptimalityTolerance',CtrlVar.Inverse.OptimalityTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FunctionTolerance',CtrlVar.Inverse.FunctionTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'StepTolerance',CtrlVar.Inverse.StepTolerance);
CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',Hfunc);

if CtrlVar.Inverse.MatlabOptimisationHessianParameters.Algorithm=="trust-region-reflective"



    switch  CtrlVar.Inverse.HessianOptions

        case "-DirectAdjoint-"
            %% Hessian provided

            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn','objective');
            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianMultiplyFcn',[]);

            Aeq = [];
            beq = [];
            nonlcon = [];


            [p,J,exitflag,output,lambda,grad,hessian] = fmincon(func,p0,Aineq,bineq,Aeq,beq,plb,pub,nonlcon,CtrlVar.Inverse.MatlabOptimisationHessianParameters);



        case "-FiniteDifferences-"

            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianFcn',[]);
            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessianMultiplyFcn',[]);




            if contains(CtrlVar.Inverse.MinimisationMethod,"-BandWidth")
                n=str2double(extractBetween(CtrlVar.Inverse.MinimisationMethod,"-BandWidth","-")) ;
            else
                n=5;
            end

            HessPattern=spdiags(ones(1,n),-(n-1)/2:(n-1)/2,numel(p0),numel(p0));  % n-diagonal

            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'HessPattern',HessPattern);
            CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FiniteDifferenceType','central');
            % v=zeros(numel(p0),1)+1e100;
            % CtrlVar.Inverse.MatlabOptimisationHessianParameters = optimoptions(CtrlVar.Inverse.MatlabOptimisationHessianParameters,'FiniteDifferenceStepSize',v);

        otherwise

            error("CaseNotFound")

    end

else

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
            error("InversionUsingMatlabOptimisationToolbox3:InvalidParameters","CtrlVar.Inverse.MinimisationMethod invalid.")

        end

        if isfield(RunInfo.Inverse,'fmincon')
            RunInfo.Inverse.fmincon=output;
        end

    else

        fprintf('Matlab Optimisation selected, but Matlab Optimisation routine not recognized.\n')
        fprintf(' Either select fminunc or fmincon. \n')
        error(' invalid input parameters ')

    end


end


% get info about the iteration, for some reason this call is correct/OJ for both fmincon and fminunc
[stop,Outs] = fminuncOutfun();


RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;RunInfo.Inverse.Iterations(end)+Outs.iteration];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Outs.fval];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.J;Outs.StepSize];
RunInfo.Inverse.R=[RunInfo.Inverse.R;Outs.fval+NaN];
RunInfo.Inverse.I=[RunInfo.Inverse.I;Outs.fval+NaN];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;Outs.GradNorm];
RunInfo.Inverse.p=Outs.p;
% If I need some further info and want to update F




end