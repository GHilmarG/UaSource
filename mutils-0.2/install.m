function install
%INSTALL mutils package. See INSTALL and RELEASE_NOTES for details.

% Copyright 2012, Marcin Krotkiewski, University of Oslo

%% Check if we are running inside milamin, or independently
milamin_data = getappdata(0, 'milamin_data');
if ~isfield(milamin_data, 'path')
    basepath = pwd;
else
    basepath = [milamin_data.path filesep 'ext'];
end
curpath  = pwd;

%% Check if we know how to compile
try
    cd(basepath);
    fid=fopen('test.c', 'w+');
    if fid==-1
        warning([mfilename ': MATLAB MEX compiler can not be tested. Trying anyhow.']);
    else
        fprintf(fid, 'int mexFunction(){return 0;}\n');
        fclose(fid);
        mex('test.c');
        loc = dir(['test.' mexext]);
        if isempty(loc)
            warning([mfilename ': MATLAB MEX compiler is not configured correctly.'...
                ' Cannot compile a test program. MILAMIN will not use external packages.']);
            return;
        end
        delete('test.c');
        delete(['test.' mexext]);
    end
catch
    wrnstate = warning('query', 'MATLAB:DELETE:FileNotFound');
    warning('off', 'MATLAB:DELETE:FileNotFound');   
    delete('test.c');
    delete(['test.' mexext]);
    warning(wrnstate);
    warning([mfilename ': MATLAB MEX compiler is not configured correctly.'...
        ' Cannot compile test program. MILAMIN will not use external packages.']);
    cd(curpath);
    return;
end

%% Compile/Install external packages
setappdata(0, 'mutils_paths', []);
update_path([basepath]);

% triangle mesher
update_path([basepath filesep 'triangle']);
display('Installing triangle...');
try
    if triangle_install([basepath filesep 'triangle'])
        display('triangle MEX file available.');
    end
catch
    display('triangle NOT available.');
end
display(' ');

% mutils
update_path([basepath filesep 'mutils']);
display('Installing mutils...');
try
    mutils_install([basepath filesep 'mutils']);
catch
    display('mutils NOT available.');
end

% SuiteSparse
update_path([basepath filesep 'SuiteSparse']);
display('Installing needed SuiteSparse components...');
try
    if suitesparse_install([basepath filesep 'SuiteSparse'])
        display('SuiteSparse components available.');
    end
catch
    display('mutils NOT available.');
end
display(' ');

% print addpath info
display(' ');
display('--------------------------------------------------------------------');
display('Paths required by mutils have been added to your MATLAB environment.');
display('You need to either save the path from menu File->Set Path,');
display('or add the following lines to your code whenever you want to use');
display('mutils:');
display(' ');
mutils_paths = getappdata(0, 'mutils_paths');
for i=1:length(mutils_paths)
    display(['    addpath(''' mutils_paths{i} ''')']);
end

end
