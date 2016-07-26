function mutils_clean(basepath)

% Copyright 2012, Marcin Krotkiewski, University of Oslo

if nargin==0
    basepath = pwd;
end

config = mutils_config(basepath);

subdirs = {'' 'libutils' 'matlab' 'sparse' 'quadtree' 'interp'};
for i=1:numel(subdirs)
    delete([basepath filesep subdirs{i} filesep '*' config.obj_extension]);
    delete([basepath filesep subdirs{i} filesep '*.' mexext]);
end

end
