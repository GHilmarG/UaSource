function [installed, MEX_OBJ] = matlab_install(basepath)

if nargin==0
    basepath = pwd;
end
curpath = pwd;
chdir(basepath);

installed = 0;

config = mutils_config([basepath filesep '..']);

MEX_OBJ_SRC = [dir([basepath filesep '*.c']) ];
MEX_OBJ_SRC = cellfun(@strcat, repmat({[basepath filesep]}, 1, length(MEX_OBJ_SRC)), {MEX_OBJ_SRC.name},...
    'UniformOutput', false);
MEX_OBJ = regexprep(MEX_OBJ_SRC, '\.c(pp)*$', config.obj_extension);


%% check, maybe we already do have what's needed in the path
installed = 1;
for i=1:numel(MEX_OBJ)
    if ~exist(MEX_OBJ{i}, 'file')
        installed = 0;
        break;
    end
end
if installed
   chdir(curpath);
   return;
end

    
%% Compile the mex files
try   
    % compile objects
    mexflags = config.mexflags;
    mexflags{end+1} = '-c';
    for i=1:length(MEX_OBJ_SRC)
        fname =  MEX_OBJ_SRC{i};
        display(['compiling ' regexprep(fname, '\\', '\\\\')]);
        mex(mexflags{:}, config.cflags, fname);
    end
    
    installed = 1;
catch
    warning([mfilename ': compilation of ' fname ' failed.']);
end

chdir(curpath);

end
