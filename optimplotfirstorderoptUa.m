function stop = optimplotfirstorderopt(~,optimValues,state,varargin)
% OPTIMPLOTFIRSTORDEROPT Plot first-order optimality at each iteration.
%
%   STOP = OPTIMPLOTFIRSTORDEROPT(X,OPTIMVALUES,STATE) plots
%   OPTIMVALUES.firstorderopt.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTFIRSTORDEROPT as the
%   plot function
%     options = optimoptions('fmincon','PlotFcn',@optimplotfirstorderopt);
%
%   Pass the options into an optimization problem to view the plot
%      fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0],[],[],options)

%   Copyright 2006-2023 The MathWorks, Inc.

% Always return a "stop" flag of false
stop = false;

% Persistent variables
persistent thePlot

switch state
    case "iter"
        if optimValues.iteration == 1
            thePlot = matlab.internal.optimfun.plotfcns.Factory.optimplotfirstorderopt(optimValues);
            thePlot.Axes.YScale='log'    ;

        elseif optimValues.iteration > 1
            thePlot.update(optimValues);
        end
end
end
