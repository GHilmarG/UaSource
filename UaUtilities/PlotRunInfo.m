function PlotRunInfo(RunInfo)
    
    %%
    figure ;
    
    
    plot(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; ylabel(' dt ' )
    
    
    xlabel(' time ' )
    
end