

function [UserVar,RunInfo,Ruvh,Kuvh,Timing]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%% uvh assembly in parallel, using spmd over sub-domains (the elements are distributed over the workers, see BuildMuaWorkers.m)
%
% (16 Oct 2026) Revised to reduce the overheads, which previously limited the speed-up:
%
% (1) F0 and F1 are kept on the workers between calls. Only the fields that change within a uvh solve (ub, vb and h, and b and s, which
%     are calculated from h) are sent in each call. All other fields are only sent again if they have changed, eg at the start of a new
%     time step, before a new uvh solve with changed fields, or after remeshing. This is detected in each call with a fingerprint (weighted
%     sums) of all numerical fields of F0 and F1 except the changing ones, of the mesh, and of the parallel pool. Whatever changes, the
%     worst case is an unnecessary resend, not stale data on the workers. Previously F0 and F1 were sent to all workers in each call,
%     which was the largest single cost (about 0.2 s for 25,000 nodes and 8 workers, ie more than the assembly on the workers itself).
%     In the same way, CtrlVar (except CtrlVar.uvhMatrixAssembly, which is switched between calls and sent in each call), UserVar and BCs1
%     are kept on the workers, and only sent again if they have changed (compared with isequaln). Sending CtrlVar, a large nested structure,
%     took about 0.13 s, although it is only about 0.3 MB.
% (2) The assembly and the summation over the workers are done in one spmd block (previously two).
% (3) If the worker data (MUA.workers) contain the sparsity pattern (see BuildMuaWorkers.m), the workers use the faster pattern-based
%     construction of the sparse matrix, as the serial assembly does.
%
% The results are identical to those of the serial assembly (uvhMatrixAssembly.m) to round-off.
%
% (16 Oct 2026) The optional fifth output, Timing, gives the times (in seconds) of the parts of the assembly: Fingerprint (detecting changes
% of the inputs), Send (sending changed inputs to the workers), SendAll (true if F0 and F1 were sent), Spmd (the spmd block), Workers (the
% longest time on a worker within the spmd block, including the summation), Fetch (transfer of the results to the client) and Total.
% It is used by uvhAssemblyScalingTest.m. Without a parallel pool, Timing.Parallel=false and only Timing.Total is set.
%
% (17 Oct 2026) The fingerprint functions were moved to SPMDFieldsFingerprint.m (shared with uvMatrixAssemblySSTREAM_SPMD.m); the key is
% unchanged.
%
%%

narginchk(8,8)
Timing=struct("Parallel",false,"Fingerprint",0,"Send",0,"SendAll",false,"Spmd",0,"Workers",NaN,"Fetch",0,"Total",0) ;
tTotal=tic ;

persistent W   % data kept on the workers between calls (Composites), with the fingerprints and copies used to detect changes

ParPool = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already

if isempty(ParPool)

    [UserVar,RunInfo,Ruvh,Kuvh]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;

else

    if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
        CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=ParPool.NumWorkers;
    end
    nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;
    MUAworkers=MUA.workers;

    % fields that change within a uvh solve, and are sent in each call
    Changing=["ub","vb","h","b","s"] ;

    % send F0 and F1 to the workers if anything else has changed since the last call
    t0=tic ; Key=SPMDFieldsFingerprint(MUA,nW,F0,strings(0),F1,Changing) ; Timing.Fingerprint=toc(t0) ; Timing.Parallel=true ;
    SendAll = isempty(W) || ~isequal(W.Key,Key) || ~isvalid(W.Pool) || W.Pool~=ParPool ;
    Timing.SendAll=SendAll ; t0=tic ;
    if SendAll
        F0all=F0 ; F1all=F1 ;
        spmd (nW)
            F0w=F0all ; F1w=F1all ; RIw=UaRunInfo ;
        end
        W=struct() ; W.F0=F0w ; W.F1=F1w ; W.RunInfo=RIw ; W.Key=Key ; W.Pool=ParPool ;
    end

    % CtrlVar (except CtrlVar.uvhMatrixAssembly), UserVar and BCs1: sent if they have changed
    CtrlVarRest=rmfield(CtrlVar,"uvhMatrixAssembly") ;
    if SendAll || ~isfield(W,"CtrlVarRest") || ~isequaln(W.CtrlVarRest,CtrlVarRest)
        spmd (nW) ; CRw=CtrlVarRest ; end
        W.CRw=CRw ; W.CtrlVarRest=CtrlVarRest ;
    end
    if SendAll || ~isfield(W,"UserVar") || ~isequaln(W.UserVar,UserVar)
        spmd (nW) ; UVw=UserVar ; end
        W.UVw=UVw ; W.UserVar=UserVar ;
    end
    if SendAll || ~isfield(W,"BCs1") || ~isequaln(W.BCs1,BCs1)
        spmd (nW) ; BCw=BCs1 ; end
        W.BCw=BCw ; W.BCs1=BCs1 ;
    end

    Timing.Send=toc(t0) ;
    V=struct() ;
    for f=Changing ; V.(f)=F1.(f) ; end
    UMA=CtrlVar.uvhMatrixAssembly ;
    UMA.CompiledKernel=false ;   % (25 Oct 2026) the compiled kernel (CtrlVar.uvhMatrixAssembly.CompiledKernel) is not used within spmd workers, as its OpenMP threads would compete with the workers
    F0w=W.F0 ; F1w=W.F1 ; RIw=W.RunInfo ; CRw=W.CRw ; UVw=W.UVw ; BCw=W.BCw ;

    t0=tic ;
    spmd (nW)
        tw0=tic ;
        F1loc=F1w ;
        fn=fieldnames(V) ;
        for i=1:numel(fn) ; F1loc.(fn{i})=V.(fn{i}) ; end
        CtrlVarLoc=CRw ; CtrlVarLoc.uvhMatrixAssembly=UMA ;
        % (31 Oct 2026) Worker 1 packs the summed results into one structure, so that the client fetches them in a single access.
        % (Summing the matrix values as a dense vector on the global pattern, with one spmdPlus for everything, was tried and was
        % slower: the dense vector carries the whole pattern, and building the sparse matrix on the client costs more than it saves.)
        [~,~,rr,kk]=uvhMatrixAssembly(UVw,RIw,CtrlVarLoc,MUAworkers,F0w,F1loc,[],BCw);   % l1 is not used in the assembly
        Xsum=spmdPlus(rr,1) ;
        if ~isempty(kk) ; Ksum=spmdPlus(kk,1) ; else ; Ksum=[] ; end
        if spmdIndex==1
            Out=struct("X",Xsum,"K",Ksum) ;
        else
            Out=[] ;
        end
        tw=toc(tw0) ;
    end
    Timing.Spmd=toc(t0) ;
    if nargout>=5 ; Timing.Workers=max([tw{:}]) ; end
    t0=tic ;
    O=Out{1} ;   % (31 Oct 2026) a single access to the Composite (each access costs about 5 ms)
    neq=3*MUA.Nnodes ;
    Ruvh=O.X(1:neq) ;
    Kuvh=O.K ;
    Timing.Fetch=toc(t0) ;

end

Timing.Total=toc(tTotal) ;

end
