function fig=FindOrCreateFigure(FigureName,Position)

fig=findobj(0,'name',FigureName);

if isempty(fig)
    fig=figure('name',FigureName);
    if nargin>1 && ~isempty(Position) 
        fig.Position=Position;
    end
else
    fig=figure(fig);
    Position=fig.Position;
    clf(fig)
    fig.Position=Position;
    hold off
end

end