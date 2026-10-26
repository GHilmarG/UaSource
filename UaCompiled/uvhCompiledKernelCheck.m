function [ok,RK,KernelFull,KernelRonly]=uvhCompiledKernelCheck(CtrlVar,MUA)

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
% See README.md in this folder.
%
%%

persistent Printed MexFound
if isempty(Printed) ; Printed=strings(0) ; MexFound=struct() ; end

ok=false ; RK=[] ; KernelFull="" ; KernelRonly="" ; Reason="" ;

if MUA.nod>10 || MUA.nip>37
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
