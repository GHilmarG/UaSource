function SetUaPath()
%make sure that the start directory is the top directory in matlab path
UaHomeDirectory=getenv('UaHomeDirectory');

if isempty(UaHomeDirectory)
    fprintf('The environmental variable UaHomeDirectory must contain the path to Ua home directory \n' )
    fprintf('Use the setenv command to do this. \n')
    error('Ua:SetUaPath','The environmental variable UaHomeDirectory is not set')
end

return

addpath([UaHomeDirectory,'/SuiteSparse/CHOLMOD/MATLAB'],'-begin')
%addpath([UaHomeDirectory,'/Mesh2d/Mesh2d v24'],'-begin')
addpath(genpath([UaHomeDirectory,'/mutils-0.2']),'-begin')
addpath(genpath([UaHomeDirectory,'/UaUtilities']),'-begin')
%addpath(genpath([UaHomeDirectory,'/html']),'-begin')
addpath(UaHomeDirectory,'-begin')
addpath(pwd,'-begin')

end