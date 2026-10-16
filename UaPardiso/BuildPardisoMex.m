function BuildPardisoMex()

%% Compiles pardisoMex.c (the MATLAB interface to the Intel oneMKL PARDISO solver) in the folder of this file
%
%   BuildPardisoMex
%
% (16 Oct 2026) Requires a C compiler set up for MEX (mex -setup C), eg MinGW64. No MKL headers or import libraries are needed, as the
% MKL runtime library (mkl_rt) is loaded at run time (see pardisoMex.c and UaPardisoSolve.m).
%
%%

Folder=fileparts(mfilename("fullpath")) ;
mex("-R2018a","-outdir",Folder,fullfile(Folder,"pardisoMex.c")) ;
fprintf("BuildPardisoMex: %s built.\n",fullfile(Folder,"pardisoMex."+mexext)) ;

end
