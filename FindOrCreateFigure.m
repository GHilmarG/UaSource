





function [fig,FigFound]=FindOrCreateFigure(FigureName,Position,Nx,Ny)

%%
%
% Creates a new figure if a figure with the FigureName has not already been defined,
% otherwise finds the figure with the given name and makes that figure the current
% active figure.
%
% If no Position is given, or left empty, new figure will be put on somewhere on a Nx
% times Ny grid (default: Nx=4, Ny=3)
%
% Very useful utility with a proven negative grumpiness impact.
%
%
%
% Example:
%
%   FindOrCreateFigure("TestingMapping",[100 100 1000 1000]);
%
%   fig=FindOrCreateFigure("ThisIsMyFigureAndIWantToFindItAgainLaterAndReuse")  ;
%
%
%   Note:  To keep the figure frame but to clear everything else use
%
%       clf(fig)
%
% To create a pdf of figure:
%
%    exportgraphics(fig,'fig.pdf')
%
%%

persistent nFigs FigsArray

R25=isMATLABReleaseOlderThan("R2025a") ;


NumberOfInputArguments=nargin ;

if isempty(nFigs)
    nFigs=0;
end

if isempty(FigsArray)
    FigsArray=findobj('type','figure');
end

if NumberOfInputArguments<1 || isempty(FigureName)
    error("FindOrCreateFigure:NoFigureNameGiven","FigureName is a required input")
end


if NumberOfInputArguments< 3 || isempty(Nx) || isempty(Ny)
    Nx=4;
    Ny=3;
end

screensize = get( groot, 'Screensize' ) ;

% find if figure already exists and if it is already contained in the figure array
fig=[];
for I=1:numel(FigsArray)
    if ishandle(FigsArray(I))
        if strcmp(FigsArray(I).Name,FigureName)
            fig=FigsArray(I);
            break
        end
    end
end

if isempty(fig)  % does the figure exist, but is not in figure array?
    fig=findobj('type','figure','name',FigureName);
    if ~isempty(fig)
        FigsArray=findobj('type','figure'); % OK, update figure array
    end
end


%fig=findobj(0,'name',FigureName);
%fig=findobj('type','figure','name',FigureName);


if isempty(fig)
    FigFound=false;
    fig=figure('name',FigureName,NumberTitle='off');

    %% Positions
    if NumberOfInputArguments>1 && ~isempty(Position)
        if R25
            fig.Position=Position;
        end
    else


        figWidth=screensize(3)/Nx;
        figHeight=screensize(4)/Ny;
        nx=mod(floor(nFigs/Ny),Nx) ;
        ny=mod(mod(nFigs,Nx*Ny),Ny);
        % nx=0;

        if R25
            fig.OuterPosition=[nx*figWidth ny*figHeight figWidth figHeight];
        end
        nFigs=nFigs+1;
    end
    FigsArray=findobj('type','figure'); % OK, update figure array
else
    FigFound=true;
    try
        fig=figure(fig(1));
    catch
        fig=figure;
    end

    Position=fig.Position;
    if R25
        fig.Position=Position;
    end
    hold off
end


set(gcf,'Color','white')

if nargout==0
    clear fig
end


end