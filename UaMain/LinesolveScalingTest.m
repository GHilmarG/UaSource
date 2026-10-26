function Results=LinesolveScalingTest(FileName,Options)

%% Performance of the linear solvers (MATLAB and PARDISO) as a function of the number of cores
%
%   Results=LinesolveScalingTest(FileName)
%   Results=LinesolveScalingTest(FileName,Options)
%
% (17 Oct 2026) For a given model state (restart file), this function builds the linear system solved in each Newton-Raphson iteration
% (the uvh system of the implicit transient solver, or the uv system of the diagnostic solver), and measures the time of its solution
% with MATLAB's sparse direct solver and with PARDISO (Intel oneMKL, see UaPardisoSolve.m), using different numbers of cores.
%
% Cores, not workers: the linear system is solved by the client, not by the workers of a parallel pool, which are idle meanwhile. Both
% solvers use threads: MATLAB's solver multithreaded BLAS, PARDISO the threads of the TBB layer of oneMKL. PARDISO with the TBB layer
% does not allow its number of threads to be limited from within MATLAB (MKL_Set_Num_Threads has no effect), and therefore each
% measurement is done in a separate MATLAB process that is restricted to the given number of cores by the operating system (processor
% affinity: start /affinity on Windows, taskset on Linux). Within this process both solvers see only these cores. One logical processor
% is used per physical core, filling one processor (socket) first. In addition, a measurement without restriction (all logical
% processors) and one with the sequential layer of PARDISO (one thread) are done.
%
% Measured for each number of cores:
%
%   PARDISO, each phase separately (analysis, ie reordering and symbolic factorisation; numerical factorisation; forward and back
%   substitution with iterative refinement):
%     (a) with a new analysis for each solve, ie all three phases, as MATLAB's solver always does, and
%     (b) with the analysis reused, as in Ua, where the sparsity pattern is fixed for each mesh (see UaPardisoSolve.m),
%   and the number of iterative refinement steps and of non-zeros in the factors.
%
%   MATLAB: backslash, and the factorisation (lu or chol, which include MATLAB's reordering and symbolic analysis; these can not be
%   separated) and the forward and back substitution with the factors, and the number of non-zeros in the factors.
%
%   Both solvers as called in the Newton-Raphson iterations of Ua: solveKApe.m (uvh) or solveKApeSymmetric.m (uv), including the
%   elimination of the constraints; for PARDISO with the analysis reused, and the first call (with analysis) separately.
%
% The system: for "uvh", the uvh system with its constraints (velocity boundary conditions and active thickness constraints), from
% which the constraints are eliminated as in ABfgPreEliminate.m, giving a nonsymmetric system (PARDISO matrix type 11). For "uv", the
% symmetric uv system; with nodal constraints only, the constrained unknowns are removed, giving a symmetric positive definite system
% (PARDISO matrix type 2, Cholesky), as in KKT_null_space_decomposition_chol.m; otherwise the constraints are eliminated as for uvh.
%
% A model T(p) = a + b/p + d p is fitted to the times as a function of the number of cores p (as in uvhAssemblyScalingTest.m), and a
% number of cores is recommended for PARDISO.
%
% Options (all optional):
%   .System        "uvh" (default) or "uv"
%   .Cores         numbers of cores (default 1, 2, 4, 8, ... up to the number of physical cores)
%   .Repeats       numbers of repeated measurements (the median is used): [with analysis, analysis reused] (default [3 5])
%   .MATLAB        measure MATLAB's solver as well (default true)
%   .Unrestricted  also measure without restriction of the cores (default true)
%   .Sequential    also measure PARDISO with its sequential layer (default true)
%   .PoolWorkers   open an idle thread pool with this number of workers during the measurements, to check that it does not slow down
%                  the solve (default 0, no pool)
%   .Tolerance     the recommended number of cores is the smallest number for which the time of PARDISO with the analysis reused
%                  (according to the model) is within this fraction of the smallest time (default 0.05)
%   .Plot          plot the results (default true)
%   .KeepFiles     keep the temporary files (the system, and the scripts and logs of the measurements) (default false)
%
% Results is a structure with the measurements for each case, the model parameters and the recommendation.
%
% Notes:
%   - The computer should be otherwise idle. The test takes several minutes, as each measurement is a separate MATLAB process (start-up
%     and loading of the system take 15-30 s), and MATLAB's solver takes several seconds per solve for large systems.
%   - PARDISO must be available (UA_MKL_RT, or CtrlVar.LinearSolver.PARDISO.Library in the restart file; see UaPardisoSolve.m).
%   - Windows and Linux only (macOS does not allow processor affinity to be set). On Windows, at most the first 64 logical processors
%     (one processor group) are used, and the logical processors of a physical core are assumed to be numbered consecutively; on Linux
%     the topology is read with lscpu.
%
% See the chapter "The PARDISO sparse direct solver" in the Ua Compendium.
%
% Example:
%
%   Results=LinesolveScalingTest("GreenlandExample-RestartFile.mat") ;
%   Results=LinesolveScalingTest("GreenlandExample-RestartFile.mat",struct("System","uv","Cores",[1 4 16])) ;
%
%%

if (ischar(FileName) || isstring(FileName)) && string(FileName)=="worker"   % a measurement, in a separate MATLAB process
    Worker(Options) ;
    Results=[] ;
    return
end

if nargin<2 || isempty(Options) ; Options=struct() ; end
if ismac
    error("LinesolveScalingTest:Platform","LinesolveScalingTest is only available on Windows and Linux.") ;
end

Topo=CoreTopology() ;
Opt.System="uvh" ;
Opt.Cores=unique([2.^(0:floor(log2(Topo.Physical))) Topo.Physical]) ;
Opt.Repeats=[3 5] ;
Opt.MATLAB=true ;
Opt.Unrestricted=true ;
Opt.Sequential=true ;
Opt.PoolWorkers=0 ;
Opt.Tolerance=0.05 ;
Opt.Plot=true ;
Opt.KeepFiles=false ;
fn=fieldnames(Options) ;
for i=1:numel(fn) ; Opt.(fn{i})=Options.(fn{i}) ; end
Opt.System=string(Opt.System) ;
Opt.Cores=sort(Opt.Cores(:)') ;
Opt.Cores=Opt.Cores(Opt.Cores>=1 & Opt.Cores<=numel(Topo.CoreCPUs)) ;

%% the system, built once in this session and saved for the measurements
wd=tempname ; mkdir(wd) ;
Cleanup=onCleanup(@() RemoveFolder(wd,Opt.KeepFiles)) ;
fprintf("\n LinesolveScalingTest: %s, %s system\n",FileName,Opt.System) ;
t0=tic ;
Sys=BuildSystem(FileName,Opt.System) ;
Lib=PardisoLibrary(Sys.CtrlVar) ;
SystemFile=fullfile(wd,"System.mat") ;
save(SystemFile,"Sys","-v7.3") ;
fprintf(" %s system built and saved in %.0f s: %i unknowns, %i non-zeros, %i constraints (matrix solved: %i unknowns, %i non-zeros, PARDISO matrix type %i)\n",...
    Opt.System,toc(t0),Sys.Info.n,Sys.Info.nnz,Sys.Info.m,size(Sys.A,1),nnz(Sys.A),Sys.mtype) ;
fprintf(" %i physical cores, %i logical processors\n",Topo.Physical,Topo.Logical) ;

%% the measurements, each in a separate MATLAB process
Cases=struct("Label",{},"Cores",{},"CPUs",{},"Layer",{},"DoMatlab",{}) ;
for p=Opt.Cores
    Cases(end+1)=struct("Label",sprintf("%i cores",p),"Cores",p,"CPUs",Topo.CoreCPUs(1:p),"Layer","TBB","DoMatlab",Opt.MATLAB) ; %#ok<AGROW>
end
if Opt.Unrestricted
    Cases(end+1)=struct("Label","unrestricted","Cores",NaN,"CPUs",[],"Layer","TBB","DoMatlab",Opt.MATLAB) ;
end
if Opt.Sequential
    Cases(end+1)=struct("Label","sequential","Cores",NaN,"CPUs",[],"Layer","SEQUENTIAL","DoMatlab",false) ;
end

MatlabExe=fullfile(matlabroot,"bin","matlab") ;
if ispc ; MatlabExe=MatlabExe+".exe" ; end
P=path ;
M=cell(numel(Cases),1) ;
fprintf("\n Measurements (one MATLAB process each):\n") ;
for k=1:numel(Cases)
    c=Cases(k) ;
    A=struct("SystemFile",SystemFile,"Layer",c.Layer,"Library",Lib,"Repeats",Opt.Repeats,"DoMatlab",c.DoMatlab,...
        "PoolWorkers",Opt.PoolWorkers,"Threads",c.Cores,"OutFile",fullfile(wd,sprintf("Result%i.mat",k)),"Path",P) ;
    ArgsFile=fullfile(wd,sprintf("Args%i.mat",k)) ; save(ArgsFile,"A") ;
    Script=fullfile(wd,sprintf("Measurement%i.m",k)) ;
    writelines(["S=load('"+ArgsFile+"') ; path(S.A.Path) ;" ; "LinesolveScalingTest(""worker"",'"+ArgsFile+"') ;"],Script) ;
    LogFile=fullfile(wd,sprintf("Measurement%i.log",k)) ;
    Cmd=sprintf('"%s" -batch "run(''%s'')" -logfile "%s"',MatlabExe,Script,LogFile) ;
    if ispc       % via start, also without restriction (a command beginning with a quoted path can be changed by cmd /c)
        Aff="" ;
        if ~isempty(c.CPUs) ; Mask=uint64(0) ; for cpu=c.CPUs ; Mask=bitset(Mask,cpu+1) ; end ; Aff=" /affinity "+dec2hex(Mask) ; end
        Cmd=sprintf('start ""%s /wait %s',Aff,Cmd) ;
    elseif ~isempty(c.CPUs)
        Cmd=sprintf('taskset -c %s %s',strjoin(string(c.CPUs),","),Cmd) ;
    end
    t0=tic ;
    [~,~]=system(Cmd) ;
    if isfile(A.OutFile)
        R=load(A.OutFile) ; M{k}=R.R ;
        fprintf("   %-13s done in %3.0f s (the process saw %i cores)\n",c.Label,toc(t0),M{k}.NumCores) ;
    else
        fprintf("   %-13s FAILED after %.0f s.",c.Label,toc(t0)) ;
        if isfile(LogFile) ; L=splitlines(string(fileread(LogFile))) ; L=L(strlength(strtrim(L))>0) ; fprintf(" End of the log:\n%s\n",strjoin(L(max(1,end-8):end),newline)) ; else ; fprintf("\n") ; end
    end
end

%% summary
Results=struct() ;
Results.Info=struct("FileName",string(FileName),"System",Opt.System,"Unknowns",Sys.Info.n,"NonZeros",Sys.Info.nnz,"Constraints",Sys.Info.m,...
    "MatrixUnknowns",size(Sys.A,1),"MatrixNonZeros",nnz(Sys.A),"PardisoMatrixType",Sys.mtype,"PhysicalCores",Topo.Physical,...
    "LogicalProcessors",Topo.Logical,"Date",datetime("now")) ;
Results.Cases=Cases ;
Results.Measurements=M ;
ok=~cellfun(@isempty,M) ;
if ~any(ok)
    fprintf("\n No measurement succeeded.\n") ;
    return
end
Get=@(f,g) cellfun(@(r) GetField(r,f,g),M) ;
Tan=Get("Pardiso","Analysis") ; Tfa=Get("Pardiso","Factorisation") ; Tso=Get("Pardiso","Solve") ; Tall=Tan+Tfa+Tso ;
Trf=Get("Pardiso","ReuseFactorisation") ; Trs=Get("Pardiso","ReuseSolve") ; Tre=Trf+Trs ;
Tmb=Get("Matlab","Backslash") ; Tmf=Get("Matlab","Factorisation") ; Tms=Get("Matlab","Solve") ;
Tkp=Get("Pardiso","SolveKApe") ; Tkp1=Get("Pardiso","SolveKApeFirst") ; Tkm=Get("Matlab","SolveKApe") ;
Results.Times=struct("Label",string({Cases.Label}),"Cores",[Cases.Cores],"PardisoAnalysis",Tan,"PardisoFactorisation",Tfa,"PardisoSolve",Tso,...
    "PardisoAllPhases",Tall,"PardisoReuseFactorisation",Trf,"PardisoReuseSolve",Trs,"PardisoReused",Tre,"MatlabBackslash",Tmb,...
    "MatlabFactorisation",Tmf,"MatlabSolve",Tms,"SolveKApePardiso",Tkp,"SolveKApePardisoFirst",Tkp1,"SolveKApeMatlab",Tkm) ;

fprintf("\n Times in seconds (PARDISO: analysis, factorisation, substitution; MATLAB: factorisation incl. reordering, substitution):\n") ;
fprintf("   %-13s | %-38s | %-17s | %-27s | %-23s\n","","PARDISO, all three phases","PARDISO, reused","MATLAB","as called in Ua") ;
fprintf("   %-13s | %7s %7s %7s = %7s | %7s %8s | %8s %8s %8s | %10s %10s\n","cores","anal.","fact.","subst.","total","fact.","total","\","fact.","subst.","MATLAB","PARDISO") ;
for k=find(ok)'
    fprintf("   %-13s | %7.3f %7.3f %7.3f = %7.3f | %7.3f %8.3f | %8.3f %8.3f %8.3f | %10.3f %10.3f\n",Cases(k).Label,Tan(k),Tfa(k),Tso(k),Tall(k),...
        Trf(k),Tre(k),Tmb(k),Tmf(k),Tms(k),Tkm(k),Tkp(k)) ;
end
k1=find(ok,1) ;
fprintf("\n   PARDISO: %i refinement steps per solve, %.1f million non-zeros in the factors | MATLAB: %.1f million non-zeros in the factors\n",...
    M{k1}.Pardiso.RefinementSteps,M{k1}.Pardiso.FactorNonZeros/1e6,GetField(M{k1},"Matlab","FactorNonZeros")/1e6) ;
fprintf("   PARDISO as called in Ua: first call (with analysis) %.2f s, then %.3f s (%s)\n",Tkp1(k1),Tkp(k1),Cases(k1).Label) ;

% scaling with the number of cores, and the model T(p)=a+b/p+d p
iC=find(ok' & ~isnan([Cases.Cores])) ;
Model=struct() ; Rec=struct("Cores",NaN,"Text","") ;
if numel(iC)>=1
    p=[Cases(iC).Cores]' ; [~,i1]=min(p) ;
    fprintf("\n Speed-up relative to %i core(s):\n",p(i1)) ;
    fprintf("   %-13s %12s %12s %12s %12s %12s\n","cores","PARDISO all","reused","factorisation","substitution","MATLAB \") ;
    for j=1:numel(iC)
        k=iC(j) ;
        fprintf("   %-13s %12.2f %12.2f %12.2f %12.2f %12.2f\n",Cases(k).Label,Tall(iC(i1))/Tall(k),Tre(iC(i1))/Tre(k),Tfa(iC(i1))/Tfa(k),Tso(iC(i1))/Tso(k),Tmb(iC(i1))/Tmb(k)) ;
    end
    if numel(iC)>=3
        A=[ones(size(p)) 1./p p] ;
        Names=["PardisoAllPhases","PardisoReused","PardisoAnalysis","PardisoFactorisation","PardisoSolve","MatlabBackslash"] ;
        Vals={Tall,Tre,Tan,Tfa,Tso,Tmb} ;
        for j=1:numel(Names)
            y=reshape(Vals{j}(iC),[],1) ;
            if all(isfinite(y)) ; Model.(Names(j))=lsqnonneg(A,y)' ; else ; Model.(Names(j))=nan(1,3) ; end
        end
        TM=@(c,m) c(1)+c(2)./m+c(3).*m ;
        nn=1:numel(Topo.CoreCPUs) ;
        Tm=TM(Model.PardisoReused,nn) ;
        Rec.Cores=find(Tm<=(1+Opt.Tolerance)*min(Tm),1) ;
        fprintf("\n Model T(p) = a + b/p + d p (a: not divided between cores, b/p: divided, d p: growing with the number of cores):\n") ;
        for j=1:numel(Names)
            fprintf("   %-22s a=%.3f s, b=%.3f s, d=%.4f s\n",Names(j),Model.(Names(j))) ;
        end
        Rec.Text=sprintf("PARDISO with the analysis reused is within %.0f%% of its fastest time with %i cores.",100*Opt.Tolerance,Rec.Cores) ;
        fprintf("\n %s\n",Rec.Text) ;
    end
end
Results.Model=Model ;

% comparison of the solvers
kU=find(ok' & string({Cases.Label})=="unrestricted",1) ; if isempty(kU) ; kU=iC(end) ; end
if ~isempty(kU) && isfinite(Tmb(kU))
    Results.Comparison=struct("Case",Cases(kU).Label,"AllPhases",Tmb(kU)/Tall(kU),"Reused",Tmb(kU)/Tre(kU),"InUa",Tkm(kU)/Tkp(kU),...
        "AnalysisShare",Tan(kU)/Tall(kU)) ;
    fprintf("\n MATLAB's solver compared with PARDISO (%s):\n",Cases(kU).Label) ;
    fprintf("   with all three phases (as MATLAB always does): PARDISO %.1f times faster\n",Tmb(kU)/Tall(kU)) ;
    fprintf("   with the analysis reused (as in Ua):           PARDISO %.1f times faster\n",Tmb(kU)/Tre(kU)) ;
    fprintf("   as called in Ua (incl. the elimination of the constraints): %.3f s instead of %.3f s (%.1f times faster)\n",Tkp(kU),Tkm(kU),Tkm(kU)/Tkp(kU)) ;
    fprintf("   The analysis takes %.0f%% of a solve with all three phases; reusing it is where most of PARDISO's advantage in Ua comes from.\n",100*Tan(kU)/Tall(kU)) ;
end
kS=find(ok' & string({Cases.Label})=="sequential",1) ;
if ~isempty(kS)
    fprintf("   PARDISO with its sequential layer (one thread): all three phases %.3f s, analysis reused %.3f s\n",Tall(kS),Tre(kS)) ;
end
fprintf("   (These results apply to this system and this computer. The computer should have been otherwise idle during the test.)\n\n") ;
Results.Recommendation=Rec ;

%% plot
if Opt.Plot && numel(iC)>=1
    fig=FindOrCreateFigure("LinesolveScalingTest") ; clf(fig) ;
    p=[Cases(iC).Cores] ;
    loglog(p,Tmb(iC),"ks-",MarkerFaceColor="k") ; hold on
    loglog(p,Tall(iC),"bo-",MarkerFaceColor="b") ;
    loglog(p,Tre(iC),"ro-",MarkerFaceColor="r") ;
    loglog(p,Tan(iC),"b:") ; loglog(p,Tfa(iC),"r--") ; loglog(p,Tso(iC),"r:") ;
    loglog(p,Tre(iC(1))*p(1)./p,"Color",[0.6 0.6 0.6]) ;
    xlabel("number of cores") ; ylabel("time (s)") ; xticks(p) ;
    title(sprintf("%s system: %i unknowns, %i cores",Opt.System,size(Sys.A,1),Topo.Physical)) ;
    legend(["MATLAB backslash","PARDISO, all three phases","PARDISO, analysis reused","PARDISO analysis","PARDISO factorisation",...
        "PARDISO substitution","ideal scaling"],Location="southwest") ;
end

end

%%
function Sys=BuildSystem(FileName,System)

w0=warning("off","all") ;     % eg the warning that decompositions can not be loaded
w=whos("-file",FileName) ; names=string({w.name}) ;
L=load(FileName,"F","MUA","BCs","l") ;
if ismember("CtrlVarInRestartFile",names) ; CtrlVar=getfield(load(FileName,"CtrlVarInRestartFile"),"CtrlVarInRestartFile") ; else ; CtrlVar=Ua2D_DefaultParameters() ; end
UserVar=[] ; if ismember("UserVarInRestartFile",names) ; UserVar=getfield(load(FileName,"UserVarInRestartFile"),"UserVarInRestartFile") ; end
warning(w0) ;
F=L.F ; BCs=L.BCs ; l=L.l ; MUA=L.MUA ;
D=Ua2D_DefaultParameters() ;
LibInFile="" ; if isfield(CtrlVar,"LinearSolver") && isfield(CtrlVar.LinearSolver,"PARDISO") && isfield(CtrlVar.LinearSolver.PARDISO,"Library") ; LibInFile=string(CtrlVar.LinearSolver.PARDISO.Library) ; end
CtrlVar.LinearSolver=D.LinearSolver ;     % the linear solver options of the current version (restart files may be older)
if LibInFile~="" ; CtrlVar.LinearSolver.PARDISO.Library=LibInFile ; end
CtrlVar.Parallel.uvAssembly.spmd.isOn=false ; CtrlVar.Parallel.uvhAssembly.spmd.isOn=false ;
if isfield(MUA,"dM") ; MUA.dM=[] ; end
MUA.workers=[] ;
MUA=UpdateMUA(CtrlVar,MUA) ;
MLC=BCs2MLC(CtrlVar,MUA,BCs) ;

if System=="uvh"
    CtrlVar.uvhMatrixAssembly.ZeroFields=false ; CtrlVar.uvhMatrixAssembly.Ronly=false ;
    [~,~,R,K]=uvhMatrixAssembly(UserVar,UaRunInfo,CtrlVar,MUA,F,F,l,BCs) ;
    l1=l ;
    if numel(l1.ubvb)~=numel(MLC.ubvbRhs) ; l1.ubvb=zeros(numel(MLC.ubvbRhs),1) ; end
    if numel(l1.h)~=numel(MLC.hRhs) ; l1.h=zeros(numel(MLC.hRhs),1) ; end
    [B,g0,luvh,~]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs,l1) ;
    x0=[F.ub;F.vb;F.h] ;
    f=-R-B'*luvh ; g=g0-B*x0 ;
    Pattern=MUA.uvhAssemblyPattern ;
    Nodal=false ;
elseif System=="uv"
    CtrlVar.uvAssembly.ZeroFields=false ; CtrlVar.uvMatrixAssembly.Ronly=false ;
    [~,R,K]=uvMatrixAssembly(UaRunInfo,CtrlVar,MUA,F,BCs) ;
    B=MLC.ubvbL ; g0=MLC.ubvbRhs ;
    x0=[F.ub;F.vb] ;
    f=-R ;
    if isempty(B) ; B=sparse(0,numel(x0)) ; g=zeros(0,1) ; else ; g=g0-B*x0 ; end
    Pattern=MUA.uvAssemblyPattern ;
    Nodal=all(full(sum(B~=0,2))==1) ;
else
    error("LinesolveScalingTest:System","Options.System must be ""uvh"" or ""uv"".") ;
end

n=size(K,1) ;
if Nodal
    % symmetric positive definite system of the unconstrained unknowns
    [ic,jc]=find(B) ; bval=full(B(sub2ind(size(B),ic,jc))) ;
    c=zeros(size(B,1),1) ; c(ic)=jc ; xc=zeros(size(B,1),1) ; xc(ic)=g(ic)./bval ;
    fr=true(n,1) ; fr(c)=false ;
    A=K(fr,fr) ; b=f(fr)-K(fr,c)*xc ;
    mtype=2 ;
else
    % constraints eliminated as in ABfgPreEliminate.m: nonsymmetric system
    fac=1/full(mean(abs(diag(K)))) ;
    BtB=B'*B ; Q=speye(n)-BtB ;
    A=Q*(fac*K)+BtB ; b=Q*(fac*f)+B'*g ;
    mtype=11 ;
end

Sys=struct("System",System,"A",A,"b",b,"mtype",mtype,"K",K,"B",B,"f",f,"g",g,"Pattern",Pattern,"CtrlVar",CtrlVar,...
    "Info",struct("n",n,"nnz",nnz(K),"m",size(B,1))) ;

end

%%
function Worker(ArgsFile)

S=load(ArgsFile) ; A=S.A ;
setenv("UA_MKL_RT",A.Library) ;
pardisoMex('load',char(A.Library),char(A.Layer)) ;
if ~isnan(A.Threads) ; maxNumCompThreads(A.Threads) ; end
if A.PoolWorkers>0 ; evalc("parpool(""Threads"",A.PoolWorkers)") ; end
S=load(A.SystemFile) ; Sys=S.Sys ;
R=struct("NumCores",feature("numcores"),"MaxNumCompThreads",maxNumCompThreads,"Layer",string(pardisoMex('layer')),...
    "MKL",string(pardisoMex('version')),"Pardiso",struct(),"Matlab",struct()) ;
fprintf("worker: %i cores, layer %s\n",R.NumCores,R.Layer) ;
Mat=Sys.A ; b=Sys.b ; mtype=Sys.mtype ;

% PARDISO, all three phases (a new analysis for each solve)
T=nan(A.Repeats(1),3) ;
for r=1:A.Repeats(1)
    h=pardisoMex('init',mtype) ;
    t0=tic ; pardisoMex('analyse',h,Mat) ; T(r,1)=toc(t0) ;
    t0=tic ; pardisoMex('factor',h,Mat) ; T(r,2)=toc(t0) ;
    t0=tic ; xP=pardisoMex('solve',h,Mat,b) ; T(r,3)=toc(t0) ;
    ip=pardisoMex('iparm',h) ;
    if r<A.Repeats(1) ; pardisoMex('free',h) ; end
end
R.Pardiso.Analysis=median(T(:,1)) ; R.Pardiso.Factorisation=median(T(:,2)) ; R.Pardiso.Solve=median(T(:,3)) ;
R.Pardiso.RefinementSteps=ip(7) ; R.Pardiso.FactorNonZeros=ip(18) ;
R.Pardiso.MemoryKB=double([ip(15) ip(16) ip(17)]) ;   % peak memory of the analysis, permanent memory, memory of the factorisation and solution (KB)

% PARDISO, analysis reused (new values, same pattern)
T=nan(A.Repeats(2),2) ;
for r=1:A.Repeats(2)
    Mr=Mat*(1+1e-12*r) ;
    t0=tic ; pardisoMex('factor',h,Mr) ; T(r,1)=toc(t0) ;
    t0=tic ; x=pardisoMex('solve',h,Mr,b) ; T(r,2)=toc(t0) ; %#ok<NASGU>
end
pardisoMex('free',h) ;
R.Pardiso.ReuseFactorisation=median(T(:,1)) ; R.Pardiso.ReuseSolve=median(T(:,2)) ;
fprintf("worker: PARDISO done\n") ;

% both solvers as called in the Newton-Raphson iterations of Ua (incl. the elimination of the constraints)
CP=Sys.CtrlVar ; CP.LinearSolver.Type="PARDISO" ; CP.LinearSolver.PARDISO.ThreadingLayer=A.Layer ;
n=size(Sys.K,1) ; m=size(Sys.B,1) ; z=zeros(n,1) ; zl=zeros(m,1) ;
if Sys.System=="uvh"
    Solve=@(C) solveKApe(Sys.K,Sys.B,Sys.f,Sys.g,z,zl,C,[],Sys.Pattern) ;
else
    Solve=@(C) solveKApeSymmetric(Sys.K,Sys.B,Sys.f,Sys.g,z,zl,C,[],Sys.Pattern) ;
end
t0=tic ; xKP=Solve(CP) ; R.Pardiso.SolveKApeFirst=toc(t0) ;
T=nan(A.Repeats(2),1) ;
for r=1:A.Repeats(2) ; t0=tic ; xKP=Solve(CP) ; T(r)=toc(t0) ; end
R.Pardiso.SolveKApe=median(T) ;

% MATLAB
R.Matlab=struct("Backslash",NaN,"Factorisation",NaN,"Solve",NaN,"FactorNonZeros",NaN,"SolveKApe",NaN,"Agreement",NaN,"AgreementInUa",NaN) ;
if A.DoMatlab
    T=nan(A.Repeats(1),1) ;
    for r=1:A.Repeats(1) ; t0=tic ; xM=Mat\b ; T(r)=toc(t0) ; end
    R.Matlab.Backslash=median(T) ;
    if mtype==2
        t0=tic ; [Rc,flag,pp]=chol(Mat,"vector") ; R.Matlab.Factorisation=toc(t0) ;
        if flag==0
            t0=tic ; xc=zeros(size(b)) ; xc(pp)=Rc\(Rc'\b(pp)) ; R.Matlab.Solve=toc(t0) ; %#ok<NASGU>
            R.Matlab.FactorNonZeros=nnz(Rc) ;
        end
        clear Rc
    else
        t0=tic ; [Lf,Uf,Pf,Qf,Rf]=lu(Mat) ; R.Matlab.Factorisation=toc(t0) ;
        t0=tic ; xl=Qf*(Uf\(Lf\(Pf*(Rf\b)))) ; R.Matlab.Solve=toc(t0) ; %#ok<NASGU>
        R.Matlab.FactorNonZeros=nnz(Lf)+nnz(Uf) ;
        clear Lf Uf Pf Qf Rf
    end
    CM=Sys.CtrlVar ; CM.LinearSolver.Type="MATLAB" ;
    T=nan(min(2,A.Repeats(1)),1) ;
    for r=1:numel(T) ; t0=tic ; xKM=Solve(CM) ; T(r)=toc(t0) ; end
    R.Matlab.SolveKApe=median(T) ;
    R.Matlab.Agreement=norm(xP-xM)/norm(xM) ;
    R.Matlab.AgreementInUa=norm(xKP-xKM)/norm(xKM) ;
    fprintf("worker: MATLAB done\n") ;
end

save(A.OutFile,"R") ;
if A.PoolWorkers>0 ; delete(gcp("nocreate")) ; end

end

%%
function T=CoreTopology()
% physical cores, logical processors, and the first logical processor of each physical core (0-based), one processor (socket) first

T.Physical=feature("numcores") ;
if ispc
    T.Logical=str2double(getenv("NUMBER_OF_PROCESSORS")) ;
    if isnan(T.Logical) ; T.Logical=T.Physical ; end
    tpc=max(1,round(T.Logical/T.Physical)) ;      % logical processors per core (assumed to be numbered consecutively)
    T.CoreCPUs=(0:T.Physical-1)*tpc ;
    T.CoreCPUs=T.CoreCPUs(T.CoreCPUs<64) ;        % affinity masks: one processor group
else
    [s,o]=system("lscpu -p=CPU,CORE,SOCKET") ;
    if s~=0 ; error("LinesolveScalingTest:lscpu","lscpu is needed to determine the cores.") ; end
    L=splitlines(string(o)) ; L=L(strlength(L)>0 & ~startsWith(L,"#")) ;
    V=double(split(L,",")) ;
    T.Logical=size(V,1) ;
    [~,i]=unique(V(:,[3 2]),"rows","first") ;      % one logical processor per (socket, core), sorted by socket and core
    T.CoreCPUs=V(i,1)' ;
    T.Physical=numel(T.CoreCPUs) ;
end

end

%%
function Lib=PardisoLibrary(CtrlVar)
Lib=string(getenv("UA_MKL_RT")) ;
if Lib=="" && isfield(CtrlVar.LinearSolver,"PARDISO") && isfield(CtrlVar.LinearSolver.PARDISO,"Library")
    Lib=string(CtrlVar.LinearSolver.PARDISO.Library) ;
end
if Lib=="" || ~isfile(Lib)
    error("LinesolveScalingTest:PARDISO","The oneMKL library was not found: set UA_MKL_RT (see UaPardisoSolve.m).") ;
end
end

%%
function v=GetField(R,Group,Name)
v=NaN ;
if isstruct(R) && isfield(R,Group) && isfield(R.(Group),Name) && ~isempty(R.(Group).(Name)) ; v=R.(Group).(Name) ; end
end

%%
function RemoveFolder(wd,Keep)
if Keep
    fprintf(" The temporary files are kept in %s\n",wd) ;
elseif isfolder(wd)
    rmdir(wd,"s") ;
end
end
