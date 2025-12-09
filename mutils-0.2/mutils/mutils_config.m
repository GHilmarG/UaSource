function config = mutils_config(basepath)

% Copyright 2012, Marcin Krotkiewski, University of Oslo

config = [];

if ispc
    config.cflags = 'COMPFLAGS=$COMPFLAGS';
    config.cxxflags = 'COMPFLAGS=$COMPFLAGS';
    config.ldflags = 'LINKFLAGS=$LINKFLAGS';
    config.obj_extension = '.obj' ;
else
    config.cflags = 'CFLAGS=\$CFLAGS';
    config.cxxflags = 'CXXFLAGS=\$CXXFLAGS';
    config.ldflags = 'LDFLAGS=\$LDFLAGS';
    config.obj_extension = '.o' ;
end

curpath = pwd;
cd(basepath);

if ~exist('mutils_compiler')
    cc = mex('mutils_compiler.c', '-output', 'mutils_compiler');
end
cc = mutils_compiler;

% compiler
config.cflags = [config.cflags ' -I"' basepath '"'];

% linker flags
if strcmp(cc, 'cl')
     %config.cflags = [config.cflags ' -DUSE_OPENMP /openmp '];
     %config.ldflags = [config.ldflags ' /openmp'];
end

if strcmp(cc, 'icc')
    config.cflags = [config.cflags ' -Wall'];
    config.cflags = [config.cflags ' -DUSE_OPENMP -openmp'];
    config.ldflags = [config.ldflags ' -openmp'];
end

if strcmp(cc, 'gcc')
    config.cflags = [config.cflags ' -Wall'];
    config.cflags = [config.cflags ' -DUSE_OPENMP -fopenmp'];
    config.ldflags = [config.ldflags ' -fopenmp '];
end

config.mexflags{1} = '-O';
if (~isempty (strfind (computer, '64')))
    config.mexflags{end+1} = '-largeArrayDims';
end
% config.mexflags{end+1} = '-v';

cd(curpath);

end
