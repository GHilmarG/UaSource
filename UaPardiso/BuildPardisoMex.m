function BuildPardisoMex()

%% Compiles pardisoMex.c (the MATLAB interface to the Intel oneMKL PARDISO solver) in the folder of this file
%
%   BuildPardisoMex
%
% (16 Oct 2026) Requires a C compiler set up for MEX (mex -setup C), eg MinGW64. No MKL headers or import libraries are needed, as the
% MKL runtime library (mkl_rt) is loaded at run time (see pardisoMex.c and UaPardisoSolve.m).
% On Linux, gcc is used, with -ldl (only needed for glibc older than 2.34).
%
%%

Folder=fileparts(mfilename("fullpath")) ;
if ismac
    error("BuildPardisoMex:Platform","The PARDISO option is not available on macOS (Intel provides no MKL for current macOS versions).")
elseif isunix
    mex("-R2018a","-outdir",Folder,fullfile(Folder,"pardisoMex.c"),"-ldl") ;   % Linux (not yet tested within MATLAB, see pardisoMex.c)
else
    mex("-R2018a","-outdir",Folder,fullfile(Folder,"pardisoMex.c")) ;
end
fprintf("BuildPardisoMex: %s built.\n",fullfile(Folder,"pardisoMex."+mexext)) ;

end
