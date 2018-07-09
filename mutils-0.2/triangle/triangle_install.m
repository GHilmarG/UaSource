function installed = triangle_install(basepath)
%TRIANGLE_INSTALL downloads and installs TRIANGLE mesh generator by Jonathan Shewchuk.
%Function returns 1 if triangle MEX function is available and working on your system.
%
%  installed = triangle_install([instpath=pwd])

% Copyright 2012, Marcin Krotkiewski, University of Oslo

if nargin==0
    basepath = pwd;
end

%% check, maybe we already do have what's needed in the path
if exist(['triangle.' mexext]) == 3
    display('triangle MEX file is already installed on this system.');
    display('Using existing MEX file.');
    installed = 1;
    return;
end

config = compiler_flags;

triangle_url = 'http://www.netlib.org/voronoi/triangle.zip';
instpath     = [basepath filesep 'sources' filesep];
curpath      = pwd;
installed    = 0;
cd(basepath);

%% download and unzip triangle sources
if ~exist(instpath)
    display([mfilename ': attempting to download and install triangle mesh generator']);
    if ~exist(instpath)
        mkdir(instpath);
    end
    ofname = [instpath 'triangle.zip'];
    [f,status] = urlwrite(triangle_url, ofname);
    if status==0
        warning([mfilename ': could not download triangle.']);
        cd(curpath);
        return;
    end
    unzip(ofname, instpath);
    delete_sources = 1;
else
    delete_sources = 0;
end

%% Compile triangle mex file
try
    cflags = [config.cflags ' -DTRILIBRARY -DNO_TIMER -I"' instpath '"'];
    mex(cflags, 'triangle_mex.c', '-output', 'triangle', [instpath 'triangle.c']);
    installed = 1;
catch
    warning([mfilename ': compilation of triangle failed.']);
end

if delete_sources
    rmdir(instpath, 's');
end
cd(curpath);

end
