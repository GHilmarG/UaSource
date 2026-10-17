

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
    t0=tic ; Key=StaticFingerprint(MUA,F0,F1,Changing,nW) ; Timing.Fingerprint=toc(t0) ; Timing.Parallel=true ;
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
    F0w=W.F0 ; F1w=W.F1 ; RIw=W.RunInfo ; CRw=W.CRw ; UVw=W.UVw ; BCw=W.BCw ;

    t0=tic ;
    spmd (nW)
        tw0=tic ;
        F1loc=F1w ;
        fn=fieldnames(V) ;
        for i=1:numel(fn) ; F1loc.(fn{i})=V.(fn{i}) ; end
        CtrlVarLoc=CRw ; CtrlVarLoc.uvhMatrixAssembly=UMA ;
        [~,~,rr,kk]=uvhMatrixAssembly(UVw,RIw,CtrlVarLoc,MUAworkers,F0w,F1loc,[],BCw);   % l1 is not used in the assembly
        rrsum = spmdPlus(rr,1);
        kksum = spmdPlus(kk,1);
        tw=toc(tw0) ;
    end
    Timing.Spmd=toc(t0) ;
    if nargout>=5 ; Timing.Workers=max([tw{:}]) ; end
    t0=tic ;

    % rrsum and kksum are composites
    % Ruvh and Kuvh are double sparse

    Ruvh=rrsum{1};
    Kuvh=kksum{1};
    Timing.Fetch=toc(t0) ;

end

Timing.Total=toc(tTotal) ;

end

%%
function Key=StaticFingerprint(MUA,F0,F1,Changing,nW)
% Fingerprint of the mesh, the number of workers, and all numerical fields of F0 and F1 except those that change within a uvh solve
% (17 Oct 2026) the key of the mesh (see MeshKey.m) instead of sum(MUA.connectivity(:)), which did not detect eg a renumbering of the
% nodes or the elements, or a change of the coordinates
if isfield(MUA,"MeshKey") && isstruct(MUA.MeshKey) && isfield(MUA.MeshKey,"Topology") && isfield(MUA.MeshKey,"Geometry")
    MK=[MUA.MeshKey.Topology MUA.MeshKey.Geometry] ;
else
    MK=MeshKey(MUA) ; MK=[MK.Topology MK.Geometry] ;
end
Key=[MUA.Nnodes MUA.Nele MK nW FieldSums(F0,strings(0)) FieldSums(F1,Changing)] ;
end

function s=FieldSums(F,Exclude)
p=setdiff(string(properties(F)),Exclude,"stable") ;
s=zeros(1,2*numel(p)) ;
for k=1:numel(p)
    x=F.(p(k)) ;
    if isstruct(x)          % eg LSFMask, GF: weighted sum over its numerical fields
        fn=fieldnames(x) ;
        for i=1:numel(fn)
            y=x.(fn{i}) ;
            if (isnumeric(y) || islogical(y)) && ~isempty(y)
                [a,b]=WeightedSum(y) ; s(2*k-1)=s(2*k-1)+i*a ; s(2*k)=s(2*k)+b ;
            end
        end
    elseif (isnumeric(x) || islogical(x)) && ~isempty(x)
        [s(2*k-1),s(2*k)]=WeightedSum(x) ;
    end
end
end

function [v,n]=WeightedSum(x)
% weighted sum of the finite values, and the number of non-finite values
x=double(full(x(:))) ; nf=~isfinite(x) ; n=nnz(nf) ; x(nf)=0 ;
w=1+mod((1:numel(x))',97)/97 ;
v=sum(x.*w) ;
end