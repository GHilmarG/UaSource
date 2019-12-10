function fig=FindOrCreateFigure(FigureName,Position)

%%
%
% Creates a new figure if a figure with the FigureName has not already been
% defined, otherwise finds the figure with the given name and makes that figure
% the current active figure.
%
% Very usefull utility with a proven negative grumpiness impact.
%
% Example:    
%
%   FindOrCreateFigure("TestingMapping",[100 100 1000 1000]);
% 
%
%%

if nargin<1 || isempty(FigureName)
    error("FindOrCreateFigure:NoFigureNameGiven","FigureName is a required input")
end


fig=findobj(0,'name',FigureName);

if isempty(fig)
    fig=figure('name',FigureName);
    if nargin>1 && ~isempty(Position) 
        fig.Position=Position;
    end
else
    fig=figure(fig);
    Position=fig.Position;
    % clf(fig)
    fig.Position=Position;
    hold off
end

end