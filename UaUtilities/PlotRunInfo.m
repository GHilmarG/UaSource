function PlotRunInfo(RunInfo)
    
    %%
    figure ;
    semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; ylabel(' dt ' ) ; xlabel(' time ' )
    title('time step as a function of time')
    
    figure ; 
    histogram(RunInfo.Forward.dt) ; xlabel('dt')
    title('dt Histogram')
    
    
    figure ; 
    plot(RunInfo.Forward.time,RunInfo.Forward.uvhIterations) ; 
    xlabel('time') ; ylabel('uvh iterations')
    title('uvh iterations as a function of time')
    
end