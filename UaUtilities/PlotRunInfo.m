function PlotRunInfo(RunInfo)
    
    %%
    
    
    FindOrCreateFigure("RunInfo uvh: time step and iterations")
    yyaxis left
    semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; 
    ylabel('time step')
    

    
    yyaxis right 
    stairs(RunInfo.Forward.time,RunInfo.Forward.uvhIterations) ; 
    ylabel('uvh iterations')

    
    xlabel('time') ; 
    legend("time step","#uvh iterations")
    
     FindOrCreateFigure("RunInfo uvh: time step histogram and iterations")
     histogram(RunInfo.Forward.dt) ; xlabel('dt')
     title('dt Histogram')
    
    
    
end