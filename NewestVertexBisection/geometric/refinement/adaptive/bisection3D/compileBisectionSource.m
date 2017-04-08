[pathstr, name, ext] = fileparts(mfilename('fullpath'));
cd(pathstr);
mex CFLAGS='$CFLAGS -Wall -std=c99 -fPIC' -O ./src/bisectionRefine3DRecursiveC.c
disp('Compilation finished');