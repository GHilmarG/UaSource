function [ok,RK,KernelFull,KernelRonly]=uvhCompiledKernelCheck(CtrlVar,MUA,Failure)

%% Decides if uvhMatrixAssembly.m can use the compiled element-based kernel (CtrlVar.uvhMatrixAssembly.CompiledKernel=true)
%
%   [ok,RK,KernelFull,KernelRonly]=uvhCompiledKernelCheck(CtrlVar,MUA)
%
% ok          : true if the compiled kernel can be used for these options and this element type
% RK          : run-time parameters for the kernel (see MakeKernelParamsFlex.m)
% KernelFull  : name of the MEX file for the residual and the Jacobian
% KernelRonly : name of the MEX file for the residual only
%
% Requires: element type with at most 10 nodes and 37 integration points, the cached assembly pattern MUA.uvhAssemblyPattern, the MEX
% files (kernel and SparseFromCSC), and options supported by the kernel (MakeKernelParamsFlex.m).
% If ok=false, uvhMatrixAssembly.m uses the vectorised assembly. The reason is printed once per MATLAB session.
%
% (26 Oct 2026) Also checks, once per session and variant, that the MEX files are up to date: built with the MATLAB release that is
% running, and none of their source files changed since (UaCompiledBuildInfo.m). Called with a third input, the error (MException)
% with which a MEX file failed, it reports the failure and switches the compiled kernel off for the rest of the session:
%
%   uvhCompiledKernelCheck(CtrlVar,MUA,ME)
%
% To check again (eg after rebuilding the MEX files): clear uvhCompiledKernelCheck
%
% See README.md in this folder.
%
%%

persistent Printed MexFound UpToDate Failed
if isempty(Printed) ; Printed=strings(0) ; MexFound=struct() ; UpToDate=struct() ; Failed="" ; end

if nargin>2 && ~isempty(Failure)   % a MEX file failed: switch the compiled kernel off for the rest of the session
    Failed="the compiled kernel failed ("+string(Failure.message)+")" ;
    fprintf(CtrlVar.fidlog,"uvhMatrixAssembly: %s. The vectorised assembly is used instead for the rest of this MATLAB session (to try again: clear uvhCompiledKernelCheck).\n",Failed) ;
    Printed(end+1)=Failed ;
    ok=false ; RK=[] ; KernelFull="" ; KernelRonly="" ;
    return
end

ok=false ; RK=[] ; KernelFull="" ; KernelRonly="" ; Reason="" ;

if Failed~=""
    Reason=Failed ;
elseif MUA.nod>10 || MUA.nip>37
    Reason=sprintf("elements with %i nodes and %i integration points are outside the compiled range (at most 10 nodes and 37 integration points)",MUA.nod,MUA.nip) ;
else
    if MUA.nod==3 && MUA.nip==6
        Variant="n3i6" ; KernelFull="uvhKernelFull_n3i6_mex" ; KernelRonly="uvhKernelRonly_n3i6_mex" ;   % fixed sizes, faster
    else
        Variant="general" ; KernelFull="uvhKernelFull_mex" ; KernelRonly="uvhKernelRonly_mex" ;
    end
    if ~isfield(MexFound,Variant)
        MexFound.(Variant) = exist(KernelFull,"file")==3 && exist(KernelRonly,"file")==3 && exist("SparseFromCSC","file")==3 ;
    end
    HasPattern = isfield(MUA,"uvhAssemblyPattern") && ~isempty(MUA.uvhAssemblyPattern) && all(isfield(MUA.uvhAssemblyPattern,["map","i0","j0","iR","TopologyKey"])) ;
    if ~HasPattern
        Reason="the mesh has no cached uvh assembly pattern (MUA.uvhAssemblyPattern, see AssemblyPatternCache.m)" ;
    elseif ~MexFound.(Variant)
        Reason="not all of the MEX files "+KernelFull+", "+KernelRonly+" and SparseFromCSC were found (they are built with BuildUaCompiled.m in UaSource/UaCompiled)" ;
    else
        if ~isfield(UpToDate,Variant)   % once per session and variant: built with this release, sources unchanged since?
            [u,Why,Note]=UaCompiledBuildInfo("check",[KernelFull,KernelRonly,"SparseFromCSC"]) ;
            UpToDate.(Variant)=struct("ok",u,"Reason",Why) ;
            if Note~="" && ~any(Printed==Note)
                Printed(end+1)=Note ;
                fprintf(CtrlVar.fidlog,"uvhMatrixAssembly: note: %s.\n",Note) ;
            end
        end
    end
    if Reason=="" && MexFound.(Variant) && ~UpToDate.(Variant).ok
        Reason=UpToDate.(Variant).Reason ;
    elseif Reason=="" && MexFound.(Variant)
        [RK,okP,Missing]=MakeKernelParamsFlex(CtrlVar) ;
        if okP
            ok=true ;
        else
            Reason="these options are not supported by the compiled kernel: "+strjoin(Missing,", ") ;
        end
    end
end

if ~ok && ~any(Printed==Reason)
    Printed(end+1)=Reason ;
    fprintf(CtrlVar.fidlog,"uvhMatrixAssembly: compiled kernel not used, as %s. The vectorised assembly is used instead.\n",Reason) ;
end

end
