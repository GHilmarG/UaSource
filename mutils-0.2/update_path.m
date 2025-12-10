function update_path(newpath)
	mutils_paths = getappdata(0, 'mutils_paths');
    mutils_paths{end+1} = newpath;
    setappdata(0, 'mutils_paths', mutils_paths);
    addpath(newpath);
end
