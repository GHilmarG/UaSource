function config = compiler_flags

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
end
