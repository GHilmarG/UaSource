function PlotRunInfo(RunInfo,FigName)

%%

if nargin< 2
    FigName="";
end

fig=FindOrCreateFigure("RunInfo: time steps and iterations"+FigName) ; clf(fig) ;
yyaxis left
semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-',DisplayName="time step") ;
ylabel('time step')

yyaxis right

I=~isnan(RunInfo.Forward.uvhIterations) ;
if numel(find(I)) > 0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.uvhIterations(I),DisplayName="#uvh iterations") ;
end

hold on

I=~isnan(RunInfo.Forward.uvIterations) ;
if numel(find(I))>0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.uvIterations(I),DisplayName="#uv iterations") ;
end

I=~isnan(RunInfo.Forward.uvIterations) ;
if numel(find(I))>0
    stairs(RunInfo.Forward.time(I),RunInfo.Forward.hIterations(I),DisplayName="#h iterations") ;
end
ylabel('iterations')

xlabel('time') ;
legend;

tt=axis; axis([tt(1) tt(2) 0 tt(4)])




FindOrCreateFigure("RunInfo: time step histogram and iterations"+FigName)
histogram(RunInfo.Forward.dt) ; xlabel('dt')
title('dt Histogram')




end