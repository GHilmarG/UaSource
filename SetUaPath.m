
function SetUaPath()

% If here, then Ua is on the path, or user is running Ua from the Ua
% installation directory.
%
% Check if the path includes, for example, Mesh2d.
% And if not, add Ua directory and all subdirectory to the path
%
if ispc
    onPath = contains(lower(path),lower('Mesh2d'));
else
    onPath = contains(path,'Mesh2d');
end

% If Mesh2d not on path, set path to include all Ua subdirectory
if ~onPath
    [filepath,name,ext]=fileparts(which('Ua'));
    addpath(genpath(filepath))
end

    
end