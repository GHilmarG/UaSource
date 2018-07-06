function [xGL,yGL]=ReadBindschadlerGroundingLine(directory)

% reads grounding line as determined by Bob
% returns coordinates in polar stereographic 

if nargin==0
    
    locdir=pwd;
    AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');
    
    if isempty(AntarcticGlobalDataSets)
        error('The environmental variable AntarcticDataSets not defined' )
    end
    
    
    cd(AntarcticGlobalDataSets);
    cd GroundingLine
else
    cd(directory)
end


load('GL','xGL','yGL')  ;
cd(locdir)

end

