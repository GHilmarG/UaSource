

%%
clear all ; close all

TestFunctionParameters.Name=string('sphere') ;


CtrlVar.Inverse.InitialLineSearchStepSize=1;
CtrlVar.Inverse.HessianEstimate='0';
CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize=1e-5;
CtrlVar.Inverse.MinimumRelativelLineSearchStepSize=1e-4;
CtrlVar.Inverse.MaximumNumberOfLineSeachSteps=100;
CtrlVar.Inverse.InfoLevel=1;
CtrlVar.Inverse.Iterations=20;
CtrlVar.Inverse.GradientUpgradeMethod='conjgrad' ; %{'SteepestDecent','ConjGrad'}
CtrlVar.ConjugatedGradientsUpdate='FR';
%CtrlVar.Inverse.GradientUpgradeMethod='steepestdecent' ; %{'SteepestDecent','ConjGrad'}
CtrlVar.NewtonAcceptRatio=0.5;
CtrlVar.NLtol=1e-15;
CtrlVar.doplots=1;
CtrlVar.Inverse.StoreSolutionAtEachIteration=1;
CtrlVar.ConjugatedGradientsRestartThreshold=30; 
CtrlVar.Inverse.CalcGradI=true;
TestFunctionParameters.a=1;
TestFunctionParameters.b=5;
TestFunctionParameters.alpha=0;


TestFunctionParameters.Name=string('rosenbrock') ;  TestFunctionParameters.a=1; TestFunctionParameters.b=10;
%TestFunctionParameters.Name="ellipse" ; TestFunctionParameters.a=1; TestFunctionParameters.b=4;
TestFunctionParameters.Name=string('matyas') ; TestFunctionParameters.a=0.26; TestFunctionParameters.b=0.48;




func=@(p) TestFunction(p,TestFunctionParameters) ;

p0=[0 ; -1];
func(p0)


RunInfo=UaRunInfo;

[p,RunInfo]=UaOptimisation(CtrlVar,func,p0,RunInfo);



x=linspace(-1.5,1.5);
y=linspace(-1.5,1.5);
Z=zeros(numel(x),numel(y));

for I=1:numel(x)
    for J=1:numel(y)
        Z(I,J)=func([x(I) y(J)]);
        
    end
end


p=zeros(numel(RunInfo.Inverse.p),2);
for I=1:numel(RunInfo.Inverse.p)
    p(I,1)=RunInfo.Inverse.p{I}(1);
    p(I,2)=RunInfo.Inverse.p{I}(2);
end

Fig1=figure ; contourf(x,y,Z',40,'LineStyle','none')
colorbar
hold on
plot(p(:,1), p(:,2),'-or')

axis equal

%%


CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'FunctionTolerance',1e-10,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI,...
    'HessianFcn','objective');

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI);
% end, MatlabOptimisation parameters.




fminuncOutfun(CtrlVar);

[p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,RunInfo);


p=zeros(numel(RunInfo.Inverse.p),2);
for I=1:numel(RunInfo.Inverse.p)
    p(I,1)=RunInfo.Inverse.p{I}(1);
    p(I,2)=RunInfo.Inverse.p{I}(2);
end

figure(Fig1) 
hold on

plot(p(:,1), p(:,2),'-xg')

legend('Rosenbrock','Ua:conjgrad','Matlab:Quasi-Newton')
xlabel('x') ; ylabel('y')
