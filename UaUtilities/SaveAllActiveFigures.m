function SaveAllActiveFigures(FileName)
    
    % Finds all active figures and saves them in a file
    
    if nargin==0 || isempty(FileName)
        FileName="AllFigures.fig";
    end
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    fprintf(' Saving all active figures into the file %s ...',FileName)
    savefig(FigList,FileName)
    fprintf('...done!')
end