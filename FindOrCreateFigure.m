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
    
    persistent nFigs
    
    
    if isempty(nFigs)
        nFigs=0;
    end
    
    if nargin<1 || isempty(FigureName)
        error("FindOrCreateFigure:NoFigureNameGiven","FigureName is a required input")
    end
    
    screensize = get( groot, 'Screensize' ) ;
    
    fig=findobj(0,'name',FigureName);
    
    if isempty(fig)
        fig=figure('name',FigureName);
        if nargin>1 && ~isempty(Position)
            fig.Position=Position;
        else
            
            Nx=4 ; Ny=3;
            figWidth=screensize(3)/Nx;
            figHeight=screensize(4)/Ny;
            nx=mod(floor(nFigs/Ny),Nx) ;
            ny=mod(mod(nFigs,Nx*Ny),Ny);
            
            
            
            fig.OuterPosition=[nx*figWidth ny*figHeight figWidth figHeight];
            nFigs=nFigs+1;
        end
    else
        fig=figure(fig);
        Position=fig.Position;
        clf(fig)
        fig.Position=Position;
        hold off
    end
    
    
    
end