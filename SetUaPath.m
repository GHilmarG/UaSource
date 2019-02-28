
function SetUaPath()

%make sure that the start directory is the top directory in matlab path
UaHomeDirectory=getenv('UaHomeDirectory');

if isempty(UaHomeDirectory)
    fprintf('The environmental variable UaHomeDirectory is not set. \n' )
    fprintf('Use the setenv command to do this. \n')
    
    locdir=pwd;
    Is=strfind(locdir,filesep+"Ua"+filesep);
    if ~isempty(Is)
        cd(locdir(1:Is))
        UaFile=dir(fullfile(pwd,'**','Ua.m'));
        if ~isempty(UaFile)
            fprintf('The file Ua.m is found in the local directory. \n')
            fprintf('Will assume that this is the Ua installation directory. \n')
            fprintf('The folder %s is added to the matlab path.\n',UaFile.folder)
            UaHomeDirectory=UaFile.folder;
            setenv('UaHomeDirectory',UaHomeDirectory);
            addpath(genpath(UaHomeDirectory))
        else
            error('Ua:SetUaPath','The environmental variable UaHomeDirectory is not set')
        end
    end
    
%     
%     addpath([UaHomeDirectory,'/SuiteSparse/CHOLMOD/MATLAB'],'-begin')
%     %addpath([UaHomeDirectory,'/Mesh2d/Mesh2d v24'],'-begin')
%     addpath(genpath([UaHomeDirectory,'/mutils-0.2']),'-begin')
%     addpath(genpath([UaHomeDirectory,'/UaUtilities']),'-begin')
%     %addpath(genpath([UaHomeDirectory,'/html']),'-begin')
%     addpath(UaHomeDirectory,'-begin')
%     addpath(pwd,'-begin')
    
end