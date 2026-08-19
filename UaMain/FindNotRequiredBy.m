function unrequired = FindNotRequiredBy(targetFile, folder)
% FIND_NOT_REQUIRED_BY  List .m files in folder not required by targetFile.
%   unrequired = FIND_NOT_REQUIRED_BY(targetFile) uses current folder.
%   unrequired = FIND_NOT_REQUIRED_BY(targetFile, folder) uses specified folder.
%
%   Inputs:
%     targetFile - path to the m-file that may require others
%     folder     - (optional) folder to search for .m files (default: pwd)
%
%   Output:
%     unrequired - cell array of filenames (with .m) in folder that are NOT
%                  among the static dependencies of targetFile.

if nargin < 2 || isempty(folder)
    folder = pwd;
end

% normalize inputs
targetFile = char(targetFile);
folder = char(folder);

% list .m files in folder
d = dir(fullfile(folder,'*.m'));
files = {d.name}';

% get required files for targetFile
try
    req = matlab.codetools.requiredFilesAndProducts(targetFile);
catch
    error('Could not analyze dependencies for %s', targetFile)
end

% normalize required names to just filenames in that folder
reqNames = cellfun(@(p) lower(getFilenameIfInFolder(p, folder)), req, 'UniformOutput', false);
reqNames = reqNames(~cellfun(@isempty, reqNames)); % drop outside-folder entries
reqNames = unique(reqNames);

targetBase = lower(getFilename(targetFile));
%unrequired = files(~ismember(lower(files), [reqNames; {targetBase}]));
unrequired = files(~ismember(lower(files), reqNames));

end

function fn = getFilenameIfInFolder(p, folder)
% return filename.ext if p is inside folder, otherwise empty
p = char(p);
[pp, n, e] = fileparts(p);
fn = '';
if ~isempty(pp)
    try
        if strcmpi(fullfile(pp,''), fullfile(folder,''))
            fn = [n e];
        end
    catch
        fn = '';
    end
end
end

function fn = getFilename(p)
[~, n, e] = fileparts(p);
fn = [n e];
end
