function PlotRunInfo(RunInfo,FigName)
    
    %%
    
    if nargin< 2
        FigName="";
    end

    fig=FindOrCreateFigure("RunInfo: time steps and iterations"+FigName) ; clf(fig) ;
    yyaxis left
    semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; 
    ylabel('time step')
    
    yyaxis right 
    stairs(RunInfo.Forward.time,RunInfo.Forward.uvhIterations) ; 
    ylabel('iterations')
    hold on 
    stairs(RunInfo.Forward.time,RunInfo.Forward.uvIterations) ; 
    
    xlabel('time') ; 
    legend("time step","#uvh iterations","#uv iterations")
   
    tt=axis; axis([tt(1) tt(2) 0 tt(4)])
    

  

     FindOrCreateFigure("RunInfo: time step histogram and iterations"+FigName)
     histogram(RunInfo.Forward.dt) ; xlabel('dt')
     title('dt Histogram')
    
     
     
    
end