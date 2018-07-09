function mutils_install(basepath)

% Copyright 2012, Marcin Krotkiewski, University of Oslo
curpath=pwd;
cd(basepath);

mutils_clean(basepath);

if nargin==0
    basepath = pwd;
end

update_path([basepath filesep 'libutils']);
if libutils_install([basepath filesep 'libutils']);
    display('libutils compiled.');
end
display(' ');

update_path([basepath filesep 'matlab']);
if matlab_install([basepath filesep 'matlab']);
    display('matlab utilities compiled.');
end
display(' ');

update_path([basepath filesep 'sparse']);
if sparse_install([basepath filesep 'sparse']);
    display('sparse_create available.');
end
display(' ');

update_path([basepath filesep 'quadtree']);
if quadtree_install([basepath filesep 'quadtree']);
    display('quadtree available.');
end
display(' ');

update_path([basepath filesep 'interp']);
if interp_install([basepath filesep 'interp']);
    display('einterp available.');
end
display(' ');

delete(['mutils_compiler.' mexext]);

cd(curpath);

end
