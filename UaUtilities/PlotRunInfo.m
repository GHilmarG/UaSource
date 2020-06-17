function PlotRunInfo(RunInfo)
    
    %%
    figure ;
    semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; ylabel(' dt ' ) ; xlabel(' time ' )
    
    
    figure ; 
    histogram(RunInfo.Forward.dt) ; xlabel('dt')
    
    
    figure ; 
    plot(RunInfo.Forward.time,RunInfo.Forward.uvhIterations) ; xlabel('time') ; ylabel('uvh iterations')
    
end