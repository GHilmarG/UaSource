function [UserVar,F,l,InvFinalValues,RunInfo]=InversionUsingMatlabOptimizationToolbox2(...
    UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

%
% p is the control variable, i.e. A or C
%

switch upper(CtrlVar.Inverse.InvertFor)
    
    case 'C'
        p0=InvStartValues.C;
    case 'LOGC'
        p0=log10(InvStartValues.C);
    case {'A','LOGA'}
        p0=InvStartValues.AGlen;
end

InvFinalValues=InvStartValues;

func=@(p) JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);



%%
options = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',true);


[p,J,exitflag,output] = fminunc(func,p0,options);

[stop,Outs] = fminuncOutfun();

RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;Outs.iteration];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Outs.fval];

switch upper(CtrlVar.Inverse.InvertFor)
    
    case 'C'
        InvFinalValues.C=p;
    case 'LOGC'
        InvFinalValues.C=10.^p;
end

% If I need some further info and want to update F
%[J,Gradient,Hessian,Outs,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);


IFig1=figure('Name','Inversion','NumberTitle','off');
subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,p0) ; title('p0')
subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,p) ; title('p')
IFig1.Position=[1098.71428571429 41.5714285714286 1096 518.285714285714];


IFig2=figure('Name','True and estimated','NumberTitle','off');
subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,UserVar.TrueC) ; title('True C')
subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C) ; title('Retrieved C')
IFig2.Position=[1098.7 638.71 1096 518.29];


end