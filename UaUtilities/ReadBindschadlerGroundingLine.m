function [xGL,yGL]=ReadBindschadlerGroundingLine(directory)

% reads grounding line as determined by Bob
% returns coordinates in polar stereographic 

locdir=pwd;
if nargin==1
    cd(directory)
end


load('GroundingLinesOfAntarticaFromBobBindschadler','xGL','yGL')  ;
cd(locdir)

end

