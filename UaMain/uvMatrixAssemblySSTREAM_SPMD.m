function [Ruv,Kuv,Tint,Fext,Timing]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,BCs)

%% uv assembly in parallel, using spmd over sub-domains (the elements are distributed over the workers, see BuildMuaWorkers.m)
%
%   [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,BCs)
%   [Ruv,Kuv,Tint,Fext,Timing]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,BCs)
%
% Called by uvMatrixAssembly.m if CtrlVar.Parallel.uvAssembly.spmd.isOn=true. The results are identical to those of the serial assembly
% (uvMatrixAssemblySSTREAM.m) to round-off.
%
% (17 Oct 2026) Revised along the lines of uvhMatrixAssemblySSTREAM_SPMD.m (16 Oct 2026), to reduce the overheads:
%
% (1) F is kept on the workers between calls. Only ub and vb, which change within a uv solve, are sent in each call. All other fields
%     are only sent again if they have changed (eg AGlen or C between the forward solves of an inversion, or after remeshing), which is
%     detected with a fingerprint of the mesh, the number of workers and the other fields of F (see SPMDFieldsFingerprint.m). In the same
%     way, CtrlVar (except CtrlVar.uvAssembly, CtrlVar.uvMatrixAssembly and CtrlVar.NRuvIncomplete, which are switched between calls and
%     are sent in each call) and BCs are kept on the workers, and only sent again if they have changed (compared with isequaln).
%     Previously CtrlVar, F and BCs were sent in each call (0.17 s for 56,070 nodes and 4 thread workers, of which 0.12 s for CtrlVar).
% (2) The workers use the uv assembly pattern of their sub-domain, built once by BuildMuaWorkers.m. Previously the worker data contained
%     no uv pattern, so that each worker calculated it (and the uvh pattern) in each call (about 0.23 s in the same example).
% (3) The assembly and the summation over the workers are done in one spmd block (previously two).
% (4) Tint and Fext are summed over the workers as well, so that all outputs are as in the serial assembly (previously empty).
% (5) The number of workers is CtrlVar.Parallel.uvAssembly.spmd.nWorkers, or the size of the pool if empty, and must equal the number of
%     sub-domains in MUA.workers, which is used if they differ (previously CtrlVar.Parallel.uvhAssembly.spmd.nWorkers was set instead).
%     Without a parallel pool, or without worker data (MUA.workers), the serial assembly is used.
%
% The optional fifth output, Timing, gives the times (in seconds) of the parts of the assembly, as in uvhMatrixAssemblySSTREAM_SPMD.m:
% Fingerprint (detecting changes of the inputs), Send (sending changed inputs to the workers), SendAll (true if F was sent), Spmd (the
% spmd block), Workers (the longest time on a worker within the spmd block, including the summation), Fetch (transfer of the results to
% the client) and Total. It is used by uvAssemblyScalingTest.m. Without a parallel pool, Timing.Parallel=false and only Timing.Total is
% set.
%
%%

narginchk(4,4)
Timing=struct("Parallel",false,"Fingerprint",0,"Send",0,"SendAll",false,"Spmd",0,"Workers",NaN,"Fetch",0,"Total",0) ;
tTotal=tic ;

persistent W   % data kept on the workers between calls (Composites), with the fingerprints and copies used to detect changes

ParPool = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already

% number of sub-domains in the worker data (see BuildMuaWorkers.m); 0 if there are none or they are not valid
nBuilt=0 ;
if ~isempty(ParPool) && isfield(MUA,"workers") && ~isempty(MUA.workers)
    try
        nBuilt=nnz(exist(MUA.workers)) ;   %#ok<EXIST> exist of a Composite: which workers hold data
    catch
        nBuilt=0 ;
    end
end

if isempty(ParPool) || nBuilt==0

    [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs) ;

else

    nW=CtrlVar.Parallel.uvAssembly.spmd.nWorkers ;
    if isempty(nW) ; nW=ParPool.NumWorkers ; end
    if nW~=nBuilt ; nW=nBuilt ; end    % the workers hold the sub-domains of MUA.workers
    MUAworkers=MUA.workers ;

    % fields that change within a uv solve, and are sent in each call
    Changing=["ub","vb"] ;

    % send F to the workers if anything else has changed since the last call
    t0=tic ; Key=SPMDFieldsFingerprint(MUA,nW,F,Changing) ; Timing.Fingerprint=toc(t0) ; Timing.Parallel=true ;
    SendAll = isempty(W) || ~isequal(W.Key,Key) || ~isvalid(W.Pool) || W.Pool~=ParPool ;
    Timing.SendAll=SendAll ; t0=tic ;
    if SendAll
        Fall=F ;
        spmd (nW)
            Fw=Fall ;
        end
        W=struct() ; W.F=Fw ; W.Key=Key ; W.Pool=ParPool ;
    end

    % CtrlVar (except the fields switched between calls, which are sent in each call) and BCs: sent if they have changed
    PerCall=["uvAssembly","uvMatrixAssembly","NRuvIncomplete"] ;
    PerCall=PerCall(isfield(CtrlVar,PerCall)) ;
    CtrlVarRest=rmfield(CtrlVar,cellstr(PerCall)) ;
    PC=struct() ;
    for f=PerCall ; PC.(f)=CtrlVar.(f) ; end
    if isfield(PC,"uvMatrixAssembly") ; PC.uvMatrixAssembly.CompiledKernel=false ; end   % (30 Oct 2026) the compiled kernel is not used within spmd workers, as its OpenMP threads would compete with the workers
    if SendAll || ~isfield(W,"CtrlVarRest") || ~isequaln(W.CtrlVarRest,CtrlVarRest)
        spmd (nW) ; CRw=CtrlVarRest ; end
        W.CRw=CRw ; W.CtrlVarRest=CtrlVarRest ;
    end
    if SendAll || ~isfield(W,"BCs") || ~isequaln(W.BCs,BCs)
        spmd (nW) ; BCw=BCs ; end
        W.BCw=BCw ; W.BCs=BCs ;
    end

    Timing.Send=toc(t0) ;
    ub=F.ub ; vb=F.vb ;
    Fw=W.F ; CRw=W.CRw ; BCw=W.BCw ;

    t0=tic ;
    spmd (nW)
        tw0=tic ;
        Floc=Fw ;
        Floc.ub=ub ; Floc.vb=vb ;
        CtrlVarLoc=CRw ;
        fn=fieldnames(PC) ;
        for i=1:numel(fn) ; CtrlVarLoc.(fn{i})=PC.(fn{i}) ; end
        [rr,kk,tt,ff]=uvMatrixAssemblySSTREAM(CtrlVarLoc,MUAworkers,Floc,BCw) ;
        rrsum=spmdPlus(rr,1) ;
        kksum=spmdPlus(kk,1) ;
        ttsum=spmdPlus(tt,1) ;
        ffsum=spmdPlus(ff,1) ;
        tw=toc(tw0) ;
    end
    Timing.Spmd=toc(t0) ;
    if nargout>=5 ; Timing.Workers=max([tw{:}]) ; end
    t0=tic ;

    % rrsum, kksum, ttsum and ffsum are composites; the results are on the first worker
    Ruv=rrsum{1} ;
    Kuv=kksum{1} ;
    Tint=ttsum{1} ;
    Fext=ffsum{1} ;
    Timing.Fetch=toc(t0) ;

end

Timing.Total=toc(tTotal) ;

end
