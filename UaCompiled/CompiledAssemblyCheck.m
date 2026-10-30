function [ok,RK,KernelFull,KernelRonly]=CompiledAssemblyCheck(CtrlVar,MUA,Assembly,Failure)

%% Decides if the compiled assembly can be used, for the uv assembly (uvMatrixAssemblySSTREAM.m) or the uvh assembly (uvhMatrixAssembly.m)
%
%   [ok,RK,KernelFull,KernelRonly]=CompiledAssemblyCheck(CtrlVar,MUA,"uv")     % CtrlVar.uvMatrixAssembly.CompiledKernel=true
%   [ok,RK,KernelFull,KernelRonly]=CompiledAssemblyCheck(CtrlVar,MUA,"uvh")    % CtrlVar.uvhMatrixAssembly.CompiledKernel=true
%
% ok          : true if the compiled kernel can be used for these options and this element type
% RK          : run-time parameters for the kernel (see MakeKernelParamsFlex.m)
% KernelFull  : name of the MEX file for the residual and the Jacobian
% KernelRonly : name of the MEX file for the residual only
%
% Requires: element type with at most 10 nodes and 37 integration points, the cached assembly pattern (MUA.uvAssemblyPattern or
% MUA.uvhAssemblyPattern), the MEX files (kernel and SparseFromCSC), and options supported by the kernel (MakeKernelParamsFlex.m).
% Also checks, once per session and variant, that the MEX files are up to date: built with the MATLAB release that is running, and none
% of their source files changed since (UaCompiledBuildInfo.m).
% If ok=false, the vectorised assembly is used. The reason is printed once per MATLAB session.
%
% Called with a fourth input, the error (MException) with which a MEX file failed, it reports the failure and switches the compiled
% kernel of this assembly off for the rest of the session:
%
%   CompiledAssemblyCheck(CtrlVar,MUA,"uv",ME)
%
% uv and uvh are checked, reported and switched off separately. To check again (eg after rebuilding the MEX files):
%
%   clear CompiledAssemblyCheck
%
% See README.md in this folder.
%
% (30 Oct 2026, generalised from uvhCompiledKernelCheck.m, 25-26 Oct 2026)
%%

persistent State
if isempty(State)
    S0=struct("Printed",strings(0),"MexFound",struct(),"UpToDate",struct(),"Failed","") ;
    State=struct("uv",S0,"uvh",S0) ;
end

Assembly=string(Assembly) ;
switch Assembly
    case "uv"
        Caller="uvMatrixAssembly" ; Kernel="uvKernel" ; PatternField="uvAssemblyPattern" ; PatternFields=["map","i0","j0","iR","perm","TopologyKey"] ;
    case "uvh"
        Caller="uvhMatrixAssembly" ; Kernel="uvhKernel" ; PatternField="uvhAssemblyPattern" ; PatternFields=["map","i0","j0","iR","TopologyKey"] ;
    otherwise
        error("CompiledAssemblyCheck:Assembly","Assembly must be ""uv"" or ""uvh"".")
end
S=State.(Assembly) ;

ok=false ; RK=[] ; KernelFull="" ; KernelRonly="" ; Reason="" ;

if nargin>3 && ~isempty(Failure)   % a MEX file failed: switch the compiled kernel off for the rest of the session
    S.Failed="the compiled kernel failed ("+string(Failure.message)+")" ;
    fprintf(CtrlVar.fidlog,"%s: %s. The vectorised assembly is used instead for the rest of this MATLAB session (to try again: clear CompiledAssemblyCheck).\n",Caller,S.Failed) ;
    S.Printed(end+1)=S.Failed ;
    State.(Assembly)=S ;
    return
end

if S.Failed~=""
    Reason=S.Failed ;
elseif MUA.nod>10 || MUA.nip>37
    Reason=sprintf("elements with %i nodes and %i integration points are outside the compiled range (at most 10 nodes and 37 integration points)",MUA.nod,MUA.nip) ;
else
    if MUA.nod==3 && MUA.nip==6
        Variant="n3i6" ; KernelFull=Kernel+"Full_n3i6_mex" ; KernelRonly=Kernel+"Ronly_n3i6_mex" ;   % fixed sizes, faster
    else
        Variant="general" ; KernelFull=Kernel+"Full_mex" ; KernelRonly=Kernel+"Ronly_mex" ;
    end
    if ~isfield(S.MexFound,Variant)
        S.MexFound.(Variant) = exist(KernelFull,"file")==3 && exist(KernelRonly,"file")==3 && exist("SparseFromCSC","file")==3 ;
    end
    HasPattern = isfield(MUA,PatternField) && ~isempty(MUA.(PatternField)) && all(isfield(MUA.(PatternField),PatternFields)) ;
    if ~HasPattern
        Reason="the mesh has no cached "+Assembly+" assembly pattern (MUA."+PatternField+", see AssemblyPatternCache.m)" ;
    elseif ~S.MexFound.(Variant)
        Reason="not all of the MEX files "+KernelFull+", "+KernelRonly+" and SparseFromCSC were found (they are built with BuildUaCompiled.m in UaSource/UaCompiled)" ;
    else
        if ~isfield(S.UpToDate,Variant)   % once per session and variant: built with this release, sources unchanged since?
            [u,Why,Note]=UaCompiledBuildInfo("check",[KernelFull,KernelRonly,"SparseFromCSC"]) ;
            S.UpToDate.(Variant)=struct("ok",u,"Reason",Why) ;
            if Note~="" && ~any(S.Printed==Note)
                S.Printed(end+1)=Note ;
                fprintf(CtrlVar.fidlog,"%s: note: %s.\n",Caller,Note) ;
            end
        end
        if ~S.UpToDate.(Variant).ok
            Reason=S.UpToDate.(Variant).Reason ;
        else
            [RK,okP,Missing]=MakeKernelParamsFlex(CtrlVar) ;
            if okP
                ok=true ;
            else
                RK=[] ;
                Reason="these options are not supported by the compiled kernel: "+strjoin(Missing,", ") ;
            end
        end
    end
end

if ~ok
    KernelFull="" ; KernelRonly="" ;
    if ~any(S.Printed==Reason)
        S.Printed(end+1)=Reason ;
        fprintf(CtrlVar.fidlog,"%s: compiled kernel not used, as %s. The vectorised assembly is used instead.\n",Caller,Reason) ;
    end
end

State.(Assembly)=S ;

end
