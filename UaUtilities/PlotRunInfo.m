function [figTimeStepsAndIterations,figHist,figdtRunSteps,figIterationsRunSteps]=PlotRunInfo(RunInfo,FigName)

%%
%
% The RunInfo variable is saved in restart files and provided in DefineOutouts.
% 
% It contains information about time-step sizes, number of non-linear iterations per times step, and time-step size as a function
% of time and run steps.
%
%
% 
%
%
% Example:
%
%  PlotRunInfo(RunInfo)
%
%
%
%%

if nargin< 2
    FigName="";
end

figTimeStepsAndIterations=FindOrCreateFigure("RunInfo: time steps and iterations"+FigName) ; clf(figTimeStepsAndIterations) ;
yyaxis left
semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-',DisplayName="time step") ;
ylabel('time step, $\Delta t$',Interpreter='latex')

yyaxis right

I=~isnan(RunInfo.Forward.uvhIterations) ;
if numel(find(I)) > 0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.uvhIterations(I),DisplayName="\#uvh iterations") ;
end

hold on

I=~isnan(RunInfo.Forward.uvIterations) ;
if numel(find(I))>0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.uvIterations(I),DisplayName="\#uv iterations",LineWidth=1.2) ;
end

I=~isnan(RunInfo.Forward.hIterations) ;
if numel(find(I))>0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.hIterations(I),DisplayName="\#h iterations",LineWidth=2) ;
end


if isfield(RunInfo.Forward,"uv2hIterations")
    I=~isnan(RunInfo.Forward.uv2hIterations) ;
    if numel(find(I))>0
        stairs(RunInfo.Forward.time(I),RunInfo.Forward.uv2hIterations(I),DisplayName="\#uv-h outer iterations",LineWidth=1.5,Color="m") ;
    end
end

ylabel('\# iterations',Interpreter='latex')
xlabel('time, $t$',Interpreter='latex') ;
legend(Location="northwest",Interpreter="latex");

tt=axis; axis([tt(1) tt(2) 0 tt(4)])




figHist=FindOrCreateFigure("RunInfo: time step histogram and iterations"+FigName) ; clf(figHist) ; 

items=numel(find(~isnan( RunInfo.Forward.dt)));
nbins=max(10,fix(items/20));
histogram(RunInfo.Forward.dt,nbins,Normalization="probability") ; 
xlabel('time step, $\Delta t$',Interpreter='latex')
title("$\Delta t$ Histogram",Interpreter="latex")



figdtRunSteps=FindOrCreateFigure("RunInfo: time-steps versus run-steps"+FigName) ; clf(figdtRunSteps) ;

yyaxis left
semilogy(RunInfo.Forward.dt,'-',DisplayName="time step",LineWidth=2) ;
ylabel('time step, $\Delta t$',Interpreter='latex')

yyaxis right
plot(RunInfo.Forward.time,'-',DisplayName="time",LineWidth=2) ;
ylabel('time, $t$',Interpreter='latex')


xlabel('Run steps',Interpreter='latex')
legend(Location="northwest",Interpreter="latex");



figIterationsRunSteps=FindOrCreateFigure("RunInfo: iterations versus run-steps"+FigName) ; clf(figIterationsRunSteps) ;

yyaxis left
stairs(RunInfo.Forward.uvhIterations,'-',DisplayName="$uvh$ Iterations",LineWidth=2) ;
ylabel('$uvh$ Iterations',Interpreter='latex')

yyaxis right
stairs(RunInfo.Forward.uvIterations,'-',DisplayName="$uv$ iterations",LineWidth=2) ;
ylabel('$uv$ iterations',Interpreter='latex')


xlabel('Run Step',Interpreter='latex')
legend(Location="best",Interpreter="latex");




if ~nargout   % A trick to suppress any function output if no output requested. No need to suppress output using ;
    clearvars  figTimeStepsAndIterations figHist figdtRunSteps figIterationsRunSteps 
end





end