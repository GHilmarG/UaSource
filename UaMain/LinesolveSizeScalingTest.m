function Results=LinesolveSizeScalingTest(Files,Options)

%% Performance of the phases of PARDISO (and of MATLAB's solver) as a function of the number of unknowns
%
%   Results=LinesolveSizeScalingTest(Files)
%   Results=LinesolveSizeScalingTest(Files,Options)
%
% (20 Oct 2026) For a set of restart files of the same model with different mesh resolutions, this function measures the times of
% the three phases of PARDISO (analysis, numerical factorisation, forward and back substitution), with a new analysis for each
% solve and with the analysis reused, and of MATLAB's sparse solver, as functions of the number of unknowns n. For each file,
% LinesolveScalingTest.m is used, ie the linear system is built as in Ua, and each measurement is a separate MATLAB process
% restricted to a given number of cores.
%
% For each quantity, a power law T = c n^alpha is fitted (least squares in log-log space) to the systems with at least
% Options.MinUnknowns unknowns. For comparison, the effective exponents of n log n and n^1.5 over the same range are given: for
% two-dimensional meshes and nested-dissection orderings, the fill-in grows asymptotically as n log n, and the work of the
% factorisation as n^1.5.
%
% Files: a string array of restart files, or a file name with wildcards, eg
%
%   Results=LinesolveSizeScalingTest(fullfile(Folder,"GreenlandExample-TDR1-IR0-Stride2-nx*-RestartFile.mat")) ;
%
% Options (all optional):
%   .System        "uvh" (default) or "uv"
%   .Cores         numbers of cores of the measurements restricted by processor affinity (default 1)
%   .Unrestricted  also measure with all cores (default true)
%   .MATLAB        measure MATLAB's solver as a reference (default true)
%   .MinUnknowns   smaller systems are measured but not used in the fits (default 15000)
%   .Plot          plot the results (default true)
%
% Results: .Table (one row per file and case), .Fits (fitted exponents and prefactors per case and quantity), .Logs (the output of
% LinesolveScalingTest for each file).
%
% Notes: the computer should be otherwise idle. The test takes roughly 2 to 5 minutes per file for large systems. For systems with
% fewer than CtrlVar.LinearSolver.PARDISO.MinimumSize unknowns, UaPardisoSolve uses MATLAB's solver, so that the times "as called in
% Ua" are then those of MATLAB; the times of the phases are always measured with PARDISO directly.
%
% See also LinesolveScalingTest.m, and the chapter "The PARDISO sparse direct solver" in the Ua Compendium.
%
%%

if nargin<2 || isempty(Options) ; Options=struct() ; end
Opt.System="uvh" ;
Opt.Cores=1 ;
Opt.Unrestricted=true ;
Opt.MATLAB=true ;
Opt.MinUnknowns=15000 ;
Opt.Plot=true ;
fn=fieldnames(Options) ;
for i=1:numel(fn) ; Opt.(fn{i})=Options.(fn{i}) ; end

Files=string(Files) ;
if isscalar(Files) && contains(Files,"*")
    d=dir(Files) ;
    Files=string(fullfile({d.folder},{d.name})) ;
end
if isempty(Files) ; error("LinesolveSizeScalingTest:Files","No files given or found.") ; end
b=zeros(size(Files)) ;
for k=1:numel(Files) ; d=dir(Files(k)) ; b(k)=d.bytes ; end
[b,o]=sort(b) ; Files=Files(o) ;      % smallest first (the file size is a proxy for the number of unknowns)

Unr="" ; if Opt.Unrestricted ; Unr=" and unrestricted" ; end
fprintf("\n LinesolveSizeScalingTest: %i files, %s system, cores: %s%s, MATLAB as reference: %d\n",numel(Files),Opt.System,...
    mat2str(Opt.Cores),Unr,Opt.MATLAB) ;

Rows=struct([]) ;
Logs=strings(numel(Files),1) ;
tAll=tic ;
for k=1:numel(Files)
    % more repetitions for small systems, whose times are short
    if b(k)<10e6 ; Rep=[10 20] ; elseif b(k)<40e6 ; Rep=[5 10] ; else ; Rep=[3 5] ; end
    t0=tic ;
    R=[] ; %#ok<NASGU> assigned within evalc
    LO=struct("System",Opt.System,"Cores",Opt.Cores,"Unrestricted",Opt.Unrestricted,"Sequential",false,"MATLAB",Opt.MATLAB,...
        "Repeats",Rep,"Plot",false) ;
    try
        Logs(k)=string(evalc("R=LinesolveScalingTest(Files(k),LO) ;")) ;
    catch ME
        fprintf("   %s: FAILED (%s)\n",Files(k),ME.message) ;
        continue
    end
    [~,name]=fileparts(Files(k)) ;
    fprintf("   %-50s %8i unknowns, done in %4.0f s\n",name,R.Info.MatrixUnknowns,toc(t0)) ;
    T=R.Times ;
    for c=1:numel(R.Cases)
        M=R.Measurements{c} ;
        if isempty(M) ; continue ; end
        mem=M.Pardiso.MemoryKB ;
        Row=struct("File",string(name),"Case",string(R.Cases(c).Label),"Unknowns",R.Info.MatrixUnknowns,...
            "NonZeros",R.Info.MatrixNonZeros,"Constraints",R.Info.Constraints,"Cores",M.NumCores,...
            "FactorNonZeros",double(M.Pardiso.FactorNonZeros),"MemoryMB",max(mem(1),mem(2)+mem(3))/1024,...
            "RefinementSteps",double(M.Pardiso.RefinementSteps),...
            "Analysis",T.PardisoAnalysis(c),"Factorisation",T.PardisoFactorisation(c),"Substitution",T.PardisoSolve(c),...
            "AllPhases",T.PardisoAllPhases(c),"ReusedFactorisation",T.PardisoReuseFactorisation(c),...
            "ReusedSubstitution",T.PardisoReuseSolve(c),"Reused",T.PardisoReused(c),...
            "MatlabBackslash",T.MatlabBackslash(c),"MatlabFactorisation",T.MatlabFactorisation(c),"MatlabSubstitution",T.MatlabSolve(c),...
            "MatlabFactorNonZeros",double(GetField(M,"Matlab","FactorNonZeros")),...
            "InUaPardiso",T.SolveKApePardiso(c),"InUaMatlab",T.SolveKApeMatlab(c)) ;
        if isempty(Rows) ; Rows=Row ; else ; Rows(end+1)=Row ; end %#ok<AGROW>
    end
end
fprintf("   total time %.0f s\n",toc(tAll)) ;

Results=struct("Table",[],"Fits",[],"Logs",Logs,"Files",Files,"Options",Opt) ;
if isempty(Rows) ; return ; end
Tab=struct2table(Rows) ;
Results.Table=Tab ;

%% tables and fits, for each case
Cases=unique(Tab.Case,"stable") ;
Q=["FactorNonZeros","MemoryMB","Analysis","Factorisation","Substitution","AllPhases","ReusedFactorisation","ReusedSubstitution",...
    "Reused","MatlabBackslash"] ;
Fits=struct([]) ;
for ic=1:numel(Cases)
    S=Tab(Tab.Case==Cases(ic),:) ;
    S=sortrows(S,"Unknowns") ;
    fprintf("\n %s (the processes saw %s cores):\n",Cases(ic),strjoin(string(unique(S.Cores))',"/")) ;
    fprintf("   %9s %9s %8s | %8s %8s %8s %8s | %8s | %8s | %7s %7s\n","unknowns","fill (M)","mem (MB)","analysis","factor.",...
        "subst.","all","reused","MATLAB","M/all","M/reused") ;
    for j=1:height(S)
        fprintf("   %9i %9.2f %8.0f | %8.4f %8.4f %8.4f %8.4f | %8.4f | %8.4f | %7.2f %7.1f\n",S.Unknowns(j),S.FactorNonZeros(j)/1e6,...
            S.MemoryMB(j),S.Analysis(j),S.Factorisation(j),S.Substitution(j),S.AllPhases(j),S.Reused(j),S.MatlabBackslash(j),...
            S.MatlabBackslash(j)/S.AllPhases(j),S.MatlabBackslash(j)/S.Reused(j)) ;
    end
    use=S.Unknowns>=Opt.MinUnknowns ;
    n=S.Unknowns(use) ;
    if nnz(use)<2
        fprintf("   (fewer than two systems with at least %i unknowns: no fit)\n",Opt.MinUnknowns) ;
        continue
    end
    ln=log(n) ;
    rNlogN=polyfit(ln,log(n.*log(n)),1) ; rN15=1.5 ;
    fprintf("   Fitted exponents alpha in T = c n^alpha (%i systems, %i to %i unknowns); for comparison over this range: n log n %.2f, n^1.5 %.2f\n",...
        nnz(use),min(n),max(n),rNlogN(1),rN15) ;
    for q=Q
        y=S.(q)(use) ;
        if all(isfinite(y) & y>0)
            p=polyfit(ln,log(y),1) ;
            yf=polyval(p,ln) ; r2=1-sum((log(y)-yf).^2)/max(sum((log(y)-mean(log(y))).^2),realmin) ;
            fprintf("     %-20s alpha=%5.2f  (R^2=%.3f)\n",q,p(1),r2) ;
            F=struct("Case",Cases(ic),"Quantity",q,"Alpha",p(1),"Prefactor",exp(p(2)),"R2",r2,"MinUnknowns",min(n),"MaxUnknowns",max(n)) ;
            if isempty(Fits) ; Fits=F ; else ; Fits(end+1)=F ; end %#ok<AGROW>
        end
    end
    j=find(use,1,"last") ;
    fprintf("   Largest system (%i unknowns): the analysis is %.0f%% of a solve with all three phases; PARDISO is %.1f times faster than\n",...
        S.Unknowns(j),100*S.Analysis(j)/S.AllPhases(j),S.MatlabBackslash(j)/S.AllPhases(j)) ;
    fprintf("   MATLAB with all three phases, and %.1f times faster with the analysis reused.\n",S.MatlabBackslash(j)/S.Reused(j)) ;
end
Results.Fits=Fits ;
fprintf("\n (These results apply to these systems and this computer. The computer should have been otherwise idle during the test.)\n\n") ;

%% plot
if Opt.Plot
    fig=FindOrCreateFigure("LinesolveSizeScalingTest") ; clf(fig) ;
    tiledlayout(1,numel(Cases),"TileSpacing","compact") ;
    for ic=1:numel(Cases)
        nexttile ;
        S=sortrows(Tab(Tab.Case==Cases(ic),:),"Unknowns") ;
        n=S.Unknowns ;
        loglog(n,S.Analysis,"b-o") ; hold on
        loglog(n,S.Factorisation,"r-o") ; loglog(n,S.Substitution,"m-o") ;
        loglog(n,S.Reused,"k-s",MarkerFaceColor="k") ; loglog(n,S.MatlabBackslash,"g-^") ;
        nn=[max(Opt.MinUnknowns,min(n)) max(n)] ; y0=S.Reused(end) ;
        loglog(nn,y0*(nn/nn(2)).^1.5,"k:") ; loglog(nn,y0*(nn.*log(nn))/(nn(2)*log(nn(2))),"k--") ;
        xlabel("number of unknowns") ; ylabel("time (s)") ; title(Cases(ic)) ; grid on
        legend(["analysis","factorisation","substitution","reused (factor. + subst.)","MATLAB backslash","n^{1.5}","n log n"],Location="northwest") ;
    end
end

end

%%
function v=GetField(R,Group,Name)
v=NaN ;
if isstruct(R) && isfield(R,Group) && isfield(R.(Group),Name) && ~isempty(R.(Group).(Name)) ; v=R.(Group).(Name) ; end
end
