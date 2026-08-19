function stop = optimplotresnorm(~,optimValues,state,varargin)
    % OPTIMPLOTRESNORM Plot value of the norm of residuals at each iteration.
    %
    %   STOP = OPTIMPLOTRESNORM(X,OPTIMVALUES,STATE) plots OPTIMVALUES.resnorm.
    %
    %   Example:
    %   Create an options structure that will use OPTIMPLOTRESNORM as the plot
    %   function
    %     options = optimoptions('lsqnonlin','PlotFcn',@optimplotresnorm);
    %
    %   Pass the options into an optimization problem to view the plot
    %     lsqnonlin(@(x) sin(3*x),[1 4],[],[],options);

    %   Copyright 2006-2023 The MathWorks, Inc

    % Always return a "stop" flag of false
    stop = false;

    % Persistent variables
    persistent thePlot

    switch state
        case "iter"
            if optimValues.iteration == 0
                thePlot = matlab.internal.optimfun.plotfcns.Factory.optimplotresnorm(optimValues);
                thePlot.Axes.YScale='log'    ;

            else
                thePlot.update(optimValues);
            end

    end 
  

end
