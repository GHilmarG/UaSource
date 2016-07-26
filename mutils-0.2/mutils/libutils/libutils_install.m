function [installed LIBUTILS_OBJ] = libutils_install(basepath)

% Copyright 2012, Marcin Krotkiewski, University of Oslo

if nargin==0
    basepath = pwd;
end
curpath = pwd;
chdir(basepath);

config = mutils_config([basepath filesep '..']);

%% check, maybe we already do have what's needed in the path
% find all C files in libutils
LIBUTILS = dir([basepath filesep '*.c']);
LIBUTILS = cellfun(@strcat, repmat({[basepath filesep]}, 1, length(LIBUTILS)), {LIBUTILS.name},...
    'UniformOutput', false);
LIBUTILS_OBJ = regexprep(LIBUTILS, '\.c$', config.obj_extension);

% check if already compiled
installed = 1;
for i=1:numel(LIBUTILS_OBJ)
    if ~exist(LIBUTILS_OBJ{i}, 'file')
        installed = 0;
        break;
    end
end
if installed
   chdir(curpath);
   return;
end


%% Compile object files of LIBUTILS
try
    config.mexflags{end+1} = '-c';
    for i=1:length(LIBUTILS)
        display(['compiling ' regexprep(LIBUTILS{i}, '\\', '\\\\')]);
        mex(config.mexflags{:}, config.cflags, LIBUTILS{i});
    end
    installed = 1;
catch
    warning([mfilename ': compilation of sptools failed.']);
end

chdir(curpath);

end
