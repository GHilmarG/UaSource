function installed = sparse_install(basepath)

% Copyright 2012, Marcin Krotkiewski, University of Oslo

if nargin==0
    basepath = pwd;
end
curpath = pwd;
chdir(basepath);

config = mutils_config([basepath filesep '..']);

% list of mex functions
MEX_FUNCTIONS = {'sparse_create'};
MEX_SRC = {'sparse_create_mex.c'};
MEX_OBJ = {'opts.c'};

%% check, maybe we already do have what's needed in the path
for i=1:numel(MEX_FUNCTIONS)
    if exist([MEX_FUNCTIONS{i} '.' mexext]) == 3
        warning(['Old version of ' MEX_FUNCTIONS{i} '.' mexext ' already installed on this system will be ignored']);
    end
end


%% Compile object files of LIBUTILS
cd([basepath filesep '..' filesep 'libutils']);
[status, OBJ] = libutils_install(pwd);
cd(basepath);
MEX_OBJ = [MEX_OBJ OBJ];


%% MATLAB utility functions
cd([basepath filesep '..' filesep 'matlab']);
[status, OBJ] = matlab_install(pwd);
cd(basepath);
MEX_OBJ = [MEX_OBJ OBJ];


%% Compile the mex files
try
    for i=1:length(MEX_FUNCTIONS)
        fname =  MEX_SRC{i};
        display(['compiling ' fname]);
        ofname = MEX_FUNCTIONS{i};
        mex(config.mexflags{:}, config.cflags, fname, MEX_OBJ{:}, config.ldflags, '-output', ofname);
    end
    installed = 1;
catch
    warning([mfilename ': compilation of sptools failed.']);
end

chdir(curpath);

end
