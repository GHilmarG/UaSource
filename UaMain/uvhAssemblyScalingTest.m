function Results=uvhAssemblyScalingTest(FileName,Options)

%% Performance of the parallel (spmd) uvh assembly as a function of the number of workers
%
%   Results=uvhAssemblyScalingTest(FileName)
%   Results=uvhAssemblyScalingTest(FileName,Options)
%
% (16 Oct 2026) For a given model state, this function measures the time of the assembly of the uvh system (the residual and the
% Jacobian, and the residual only)
%
%   - with the serial assembly, using different numbers of threads, and
%   - with the parallel (spmd) assembly, using different numbers of workers n,
%
% and fits the model
%
%     T(n) = a + b/n + d n
%
% to the times of the parallel assembly. Here a is the part that does not depend on the number of workers (mainly the transfer of the
% Jacobian to the client), b/n the part that is divided between the workers (the assembly itself), and d n the part that grows with the
% number of workers (data sent to each worker, summation over the workers, competition for cores and memory bandwidth). The function
% recommends a number of workers, and estimates the gain per assembly and per Newton-Raphson iteration. For the latter, the time for the
% solution of the linear system (factorisation and solve) is also measured, and a Newton-Raphson iteration is taken to consist of one
% assembly of the residual and the Jacobian, one assembly of the residual only (line search), and one solution of the linear system.
%
% See the chapter "Parallel assembly of the uvh system: performance and scaling" in the Ua Compendium.
%
% FileName: either a file with states saved during a run with CtrlVar.uvhSaveNewtonSystem.Save=true (see SaveuvhNewtonState.m), or a
% restart file (in which case F0=F1=F).
%
% Options (all optional):
%
%   .Workers     numbers of workers (default 1, 2, 4, 8, 12, ... up to the number of cores)
%   .Threads     numbers of threads for the serial assembly (default 1, 2, 4, ... up to the number of cores)
%   .Repeats     number of repeated measurements, of which the median is used (default 5)
%   .State       index of the saved state to use (for a file with saved states; default 1)
%   .Profile     parallel profile (default "Processes")
%   .Tolerance   the recommended number of workers is the smallest number for which the time of the assembly of the residual and the
%                Jacobian (according to the model) is within this fraction of the smallest time (default 0.05)
%   .Solve       measure the solution of the linear system (default true)
%   .Plot        plot the results (default true)
%
% Results is a structure with the measurements, the model parameters and the recommendation.
%
% Notes:
%
%   - The computer should be otherwise idle during the test.
%   - The time to start the parallel pools is not included in the measurements. Starting a pool takes about 20-30 s, so the whole test
%     takes a few minutes.
%   - Parallel pools opened by this function are closed at the end. If a pool was open before, a pool with the same profile and number of
%     workers is opened again. The number of computational threads is restored.
%   - The results apply to the given mesh and computer. For small meshes the parallel assembly may be slower than the serial one.
%
% Example:
%
%   Results=uvhAssemblyScalingTest("uvhNewtonStates.mat") ;
%
%%

if nargin<2 || isempty(Options) ; Options=struct() ; end

nCores=feature('numcores') ;
Opt.Workers=unique([1 2 4 8:4:nCores nCores]) ;
Opt.Workers=Opt.Workers(Opt.Workers<=nCores) ;
Opt.Threads=unique([2.^(0:floor(log2(nCores))) nCores]) ;
Opt.Repeats=5 ;
Opt.State=1 ;
Opt.Profile="Processes" ;
Opt.Tolerance=0.05 ;
Opt.Solve=true ;
Opt.Plot=true ;
fn=fieldnames(Options) ;
for i=1:numel(fn) ; Opt.(fn{i})=Options.(fn{i}) ; end
Opt.Workers=sort(Opt.Workers(:)') ; Opt.Threads=sort(Opt.Threads(:)') ;

%% the state of the session, restored at the end (also after an error)
nThreads0=maxNumCompThreads ;
P0=gcp('nocreate') ;
if isempty(P0)
    n0=0 ; Profile0="" ;
else
    n0=P0.NumWorkers ; Profile0=string(P0.Cluster.Profile) ;
end
Cleanup=onCleanup(@() RestoreSession(nThreads0,n0,Profile0)) ;

%% the model state
[UserVar,CtrlVar,MUA,F0,F1,l1,BCs1]=LoadState(FileName,Opt.State) ;
CtrlVar.uvhMatrixAssembly.ZeroFields=false ;
CtrlVar.uvhMatrixAssembly.Ronly=false ;
CtrlVarR=CtrlVar ; CtrlVarR.uvhMatrixAssembly.Ronly=true ;
RunInfo=UaRunInfo ;

Results.Info=struct("FileName",string(FileName),"Nodes",MUA.Nnodes,"Elements",MUA.Nele,"NodesPerElement",MUA.nod,...
    "Cores",nCores,"Profile",string(Opt.Profile),"Date",datetime("now"),"Repeats",Opt.Repeats) ;

fprintf("\n uvhAssemblyScalingTest: %s\n %i nodes, %i elements (%i nodes per element), %i cores\n",FileName,MUA.Nnodes,MUA.Nele,MUA.nod,nCores) ;

%% serial assembly
fprintf("\n Serial assembly\n   threads   residual and Jacobian   residual only\n") ;
nT=numel(Opt.Threads) ;
Ser.Threads=Opt.Threads ; Ser.Full=nan(1,nT) ; Ser.ResidualOnly=nan(1,nT) ;
for k=1:nT
    maxNumCompThreads(Opt.Threads(k)) ;
    tf=nan(Opt.Repeats+1,1) ; tr=tf ;
    for r=1:Opt.Repeats+1
        F1r=Perturbed(F1,r) ;
        t0=tic ; [~,~,~,~]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1r,l1,BCs1) ; tf(r)=toc(t0) ;
        t0=tic ; [~,~,~,~]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVarR,MUA,F0,F1r,l1,BCs1) ; tr(r)=toc(t0) ;
    end
    Ser.Full(k)=median(tf(2:end)) ; Ser.ResidualOnly(k)=median(tr(2:end)) ;
    fprintf("   %5i   %12.3f s         %10.3f s\n",Opt.Threads(k),Ser.Full(k),Ser.ResidualOnly(k)) ;
end
maxNumCompThreads(nThreads0) ;
[SerialFull,iBest]=min(Ser.Full) ;
SerialResidualOnly=min(Ser.ResidualOnly) ;
Ser.BestThreads=Ser.Threads(iBest) ;
Results.Serial=Ser ;

%% solution of the linear system (factorisation and solve), with the constraints as in the uvh solver
TSolve=NaN ;
if Opt.Solve
    try
        TSolve=TimeSolve(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,Opt.Repeats) ;
        fprintf("\n Solution of the linear system (factorisation and solve, %i threads): %.3f s\n",nThreads0,TSolve) ;
    catch ME
        fprintf("\n Solution of the linear system could not be timed: %s\n",ME.message) ;
    end
end
Results.Solve=TSolve ;

%% parallel assembly
c=parcluster(Opt.Profile) ;
nW=numel(Opt.Workers) ;
Par.Workers=Opt.Workers ;
Par.Full=nan(1,nW) ; Par.ResidualOnly=nan(1,nW) ; Par.WorkerTime=nan(1,nW) ; Par.OtherSpmd=nan(1,nW) ; Par.Fetch=nan(1,nW) ;
Par.Fingerprint=nan(1,nW) ; Par.FirstCall=nan(1,nW) ; Par.MemoryPerWorkerMB=nan(1,nW) ;
wM=whos("MUA") ; wF0=whos("F0") ; wF1=whos("F1") ;

fprintf("\n Parallel assembly (profile %s)\n",Opt.Profile) ;
fprintf("   workers   residual and Jacobian (on workers, other spmd, transfer)   residual only   first call   memory per worker\n") ;
for k=1:nW
    n=Opt.Workers(k) ;
    if n>c.NumWorkers
        fprintf("   %5i   skipped: the profile %s allows at most %i workers\n",n,Opt.Profile,c.NumWorkers) ;
        continue
    end
    P=gcp('nocreate') ;
    if isempty(P) || P.NumWorkers~=n || string(P.Cluster.Profile)~=string(Opt.Profile)
        evalc("delete(P) ;") ;               % the output of delete and parpool is suppressed
        evalc("parpool(Opt.Profile,n) ;") ;
    end
    C=CtrlVar ; C.Parallel.uvhAssembly.spmd.nWorkers=n ; C.Parallel.BuildWorkers=true ;
    CR=CtrlVarR ; CR.Parallel.uvhAssembly.spmd.nWorkers=n ;
    M=MUA ; M.workers=BuildMuaWorkers(C,M,[]) ;

    % first call: all inputs are sent to the workers
    [~,~,~,~,T1]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,C,M,F0,F1,l1,BCs1) ;
    Par.FirstCall(k)=T1.Total ;

    % repeated calls, as within a uvh solve: only the changing fields are sent
    TF=repmat(T1,Opt.Repeats+1,1) ; TR=TF ;
    for r=1:Opt.Repeats+1
        F1r=Perturbed(F1,r) ;
        [~,~,~,~,TF(r)]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,C,M,F0,F1r,l1,BCs1) ;
        [~,~,~,~,TR(r)]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CR,M,F0,F1r,l1,BCs1) ;
    end
    TF=TF(2:end) ; TR=TR(2:end) ;
    Par.Full(k)=median([TF.Total]) ; Par.ResidualOnly(k)=median([TR.Total]) ;
    Par.WorkerTime(k)=median([TF.Workers]) ; Par.OtherSpmd(k)=median([TF.Spmd]-[TF.Workers]) ;
    Par.Fetch(k)=median([TF.Fetch]) ; Par.Fingerprint(k)=median([TF.Fingerprint]) ;
    Par.MemoryPerWorkerMB(k)=(wM.bytes/n+wF0.bytes+wF1.bytes)/1e6 ;
    fprintf("   %5i   %10.3f s (%.3f, %.3f, %.3f)                       %10.3f s   %8.2f s   %8.0f MB\n",n,Par.Full(k),...
        Par.WorkerTime(k),Par.OtherSpmd(k),Par.Fetch(k),Par.ResidualOnly(k),Par.FirstCall(k),Par.MemoryPerWorkerMB(k)) ;
end
evalc("delete(gcp('nocreate')) ;") ;
Results.Parallel=Par ;

%% model T(n)=a+b/n+d n, with non-negative coefficients
ok=~isnan(Par.Full) ;
n=Par.Workers(ok)' ;
Model=struct("Full",nan(1,3),"ResidualOnly",nan(1,3),"OptimalWorkers",NaN,"RmsMisfit",NaN) ;
Rec=struct("Workers",NaN,"UseParallel",false,"Full",NaN,"ResidualOnly",NaN,"SpeedUpAssembly",NaN,"GainPerNewtonIteration",NaN,"Text","") ;
if numel(n)>=3
    A=[ones(size(n)) 1./n n] ;
    Model.Full=lsqnonneg(A,Par.Full(ok)')' ;
    Model.ResidualOnly=lsqnonneg(A,Par.ResidualOnly(ok)')' ;
    Model.RmsMisfit=sqrt(mean((A*Model.Full'-Par.Full(ok)').^2)) ;
    TM=@(p,m) p(1)+p(2)./m+p(3).*m ;
    if Model.Full(3)>0
        Model.OptimalWorkers=sqrt(Model.Full(2)/Model.Full(3)) ;
    else
        Model.OptimalWorkers=Inf ;
    end
    nn=1:nCores ;
    Tm=TM(Model.Full,nn) ;
    nRec=find(Tm<=(1+Opt.Tolerance)*min(Tm),1) ;
    Rec.Workers=nRec ;
    % measured times if available, otherwise the model
    i=find(Par.Workers==nRec & ok,1) ;
    if isempty(i)
        Rec.Full=TM(Model.Full,nRec) ; Rec.ResidualOnly=TM(Model.ResidualOnly,nRec) ;
    else
        Rec.Full=Par.Full(i) ; Rec.ResidualOnly=Par.ResidualOnly(i) ;
    end
    Rec.SpeedUpAssembly=SerialFull/Rec.Full ;
    TSerialIt=SerialFull+SerialResidualOnly+TSolve ;
    TParallelIt=Rec.Full+Rec.ResidualOnly+TSolve ;
    Rec.GainPerNewtonIteration=1-TParallelIt/TSerialIt ;
    Rec.UseParallel=(Rec.Full+Rec.ResidualOnly)<(SerialFull+SerialResidualOnly) ;
else
    fprintf("\n At least three numbers of workers are needed to fit the model.\n") ;
end
Results.Model=Model ;

%% summary
fprintf("\n Model of the parallel assembly, T(n) = a + b/n + d n:\n") ;
fprintf("   residual and Jacobian: a=%.3f s, b=%.3f s, d=%.4f s (rms misfit %.3f s), optimum at n=%.1f\n",Model.Full,Model.RmsMisfit,Model.OptimalWorkers) ;
fprintf("   residual only:         a=%.3f s, b=%.3f s, d=%.4f s\n",Model.ResidualOnly) ;
fprintf("\n Best serial assembly: %.3f s (residual and Jacobian, %i threads), %.3f s (residual only)\n",SerialFull,Ser.BestThreads,SerialResidualOnly) ;
if Rec.UseParallel
    Rec.Text=sprintf("Use the parallel assembly with %i workers (CtrlVar.Parallel.uvhAssembly.spmd.isOn=true; parpool(%i)).",Rec.Workers,Rec.Workers) ;
    fprintf("\n Recommendation: %s\n",Rec.Text) ;
    fprintf("   Assembly of the residual and the Jacobian: %.3f s instead of %.3f s (%.2f times faster).\n",Rec.Full,SerialFull,Rec.SpeedUpAssembly) ;
    fprintf("   Assembly of the residual only: %.3f s instead of %.3f s.\n",Rec.ResidualOnly,SerialResidualOnly) ;
    if ~isnan(Rec.GainPerNewtonIteration)
        fprintf("   Estimated reduction of the time per Newton-Raphson iteration: %.0f%%.\n",100*Rec.GainPerNewtonIteration) ;
    end
elseif ~isnan(Rec.Workers)
    Rec.Text="Use the serial assembly: for this problem the parallel assembly is not faster on this computer." ;
    fprintf("\n Recommendation: %s\n",Rec.Text) ;
end
fprintf("   (These results apply to this mesh and this computer. The computer should have been otherwise idle during the test.)\n\n") ;
Results.Recommendation=Rec ;

%% plot
if Opt.Plot && any(ok)
    fig=FindOrCreateFigure("uvhAssemblyScalingTest") ; clf(fig) ;
    plot(Par.Workers(ok),Par.Full(ok),"ko",MarkerFaceColor="k") ; hold on
    plot(Par.Workers(ok),Par.ResidualOnly(ok),"ko") ;
    if ~any(isnan(Model.Full))
        nn=linspace(1,max(Par.Workers(ok)),200) ;
        plot(nn,Model.Full(1)+Model.Full(2)./nn+Model.Full(3).*nn,"k-") ;
        plot(nn,Model.ResidualOnly(1)+Model.ResidualOnly(2)./nn+Model.ResidualOnly(3).*nn,"k:") ;
    end
    yline(SerialFull,"--",sprintf("serial, best (%i threads)",Ser.BestThreads)) ;
    yline(SerialResidualOnly,"-.","serial, residual only") ;
    xlabel("number of workers") ; ylabel("time (s)") ; ylim([0 inf]) ;
    title(sprintf("uvh assembly: %i nodes, %i cores",MUA.Nnodes,nCores)) ;
    legend(["residual and Jacobian","residual only","model","model, residual only"],Location="northeast") ;
end

end

%%
function [UserVar,CtrlVar,MUA,F0,F1,l1,BCs1]=LoadState(FileName,iState)

w=whos("-file",FileName) ; names=string({w.name}) ;
ws=warning("off","all") ;   % eg the warning that decompositions can not be loaded
if any(startsWith(names,"State_"))
    StateNames=sort(names(startsWith(names,"State_"))) ;
    L=load(FileName,StateNames(iState)) ; S=L.(StateNames(iState)) ;
    mn=sprintf("MUA_%03d",S.MeshIndex) ;
    M=load(FileName,mn) ; MUA=M.(mn) ;
    UserVar=S.UserVar ; CtrlVar=S.CtrlVar ; F0=S.F0 ; F1=S.F1 ; l1=S.l1 ; BCs1=S.BCs1 ;
elseif all(ismember(["F","MUA","BCs"],names))
    L=struct() ; evalc("L=load(FileName,""F"",""MUA"",""BCs"") ;") ;   % evalc: the warning that decompositions (in MUA) can not be loaded is not suppressed by warning("off")
    MUA=L.MUA ; F0=L.F ; F1=L.F ; BCs1=L.BCs ;
    if ismember("l",names) ; l1=getfield(load(FileName,"l"),"l") ; else ; l1=UaLagrangeVariables ; end
    if ismember("CtrlVarInRestartFile",names)
        CtrlVar=getfield(load(FileName,"CtrlVarInRestartFile"),"CtrlVarInRestartFile") ;
    else
        CtrlVar=Ua2D_DefaultParameters() ;
    end
    if ismember("UserVarInRestartFile",names)
        UserVar=getfield(load(FileName,"UserVarInRestartFile"),"UserVarInRestartFile") ;
    else
        UserVar=struct() ;
    end
else
    warning(ws) ;
    error("uvhAssemblyScalingTest:File","%s is neither a file with saved uvh states nor a restart file.",FileName) ;
end
warning(ws) ;

if isfield(MUA,"dM") ; MUA.dM=[] ; end
MUA.workers=[] ;
if ~isfield(MUA,"uvhAssemblyPattern") || isempty(MUA.uvhAssemblyPattern)
    [~,MUA.uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA) ;
end

end

%%
function F=Perturbed(F,r)
% a small change of the velocities and the thickness, as between Newton-Raphson iterations
F.ub=F.ub*(1+1e-9*r) ; F.vb=F.vb*(1+1e-9*r) ; F.h=F.h*(1+1e-9*r) ;
end

%%
function T=TimeSolve(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,Repeats)

[~,~,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
MLC=BCs2MLC(CtrlVar,MUA,BCs1) ;
if numel(l1.ubvb)~=numel(MLC.ubvbRhs) ; l1.ubvb=zeros(numel(MLC.ubvbRhs),1) ; end
if numel(l1.h)~=numel(MLC.hRhs) ; l1.h=zeros(numel(MLC.hRhs),1) ; end
[L,cuvh,luvh,~]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1) ;
x=[F1.ub;F1.vb;F1.h] ;
if isempty(L)
    frhs=-R ; grhs=[] ;
else
    frhs=-R-L'*luvh ; grhs=cuvh-L*x ;
end
t=nan(Repeats,1) ;
for r=1:Repeats
    t0=tic ; [~,~]=solveKApe(K,L,frhs,grhs,x*0,luvh*0,CtrlVar) ; t(r)=toc(t0) ;
end
T=median(t) ;

end

%%
function RestoreSession(nThreads0,n0,Profile0)
maxNumCompThreads(nThreads0) ;
P=gcp('nocreate') ;
if n0>0
    if isempty(P) || P.NumWorkers~=n0 || string(P.Cluster.Profile)~=Profile0
        evalc("delete(P) ;") ;
        parpool(Profile0,n0) ;
    end
else
    evalc("delete(P) ;") ;
end
end