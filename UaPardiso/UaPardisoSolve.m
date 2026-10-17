function x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType,StructuralPattern)

%% Solves A x = b with the Intel oneMKL PARDISO sparse direct solver (through pardisoMex)
%
%   x=UaPardisoSolve(CtrlVar,A,b)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType,StructuralPattern)
%   Stats=UaPardisoSolve("stats")     % numbers of analyses, factorisations, solves, pattern extensions, repeated solves, fallbacks, ...
%   UaPardisoSolve("reset")           % release the memory used by PARDISO and reset the statistics
%
% (16 Oct 2026) Used by solveKApe.m, ABfgPreEliminate.m, solveKApeSymmetric.m and KKT_null_space_decomposition_chol.m if
% CtrlVar.LinearSolver.Type="PARDISO" (see Ua2D_DefaultParameters.m).
%
% MatrixType:
%   "nonsymmetric" (default)  real nonsymmetric matrix (PARDISO matrix type 11, LU factorisation)
%   "spd"                     real symmetric positive definite matrix (type 2, Cholesky factorisation). If the factorisation fails
%                             (the matrix is not positive definite), the system is solved as "symmetric" instead, with a warning.
%   "symmetric"               real symmetric indefinite matrix (type -2, LDL' factorisation with pivoting)
% For the symmetric types only the lower triangle of A is used, which in the compressed sparse column format of MATLAB is the upper
% triangle in the compressed sparse row format expected by PARDISO.
%
% Systems with fewer than CtrlVar.LinearSolver.PARDISO.MinimumSize unknowns are solved with mldivide, as for them the analysis by PARDISO
% can cost more than it saves, and diagonal systems are solved directly.
%
% Fixed pattern: The analysis (reordering and symbolic factorisation) is much more expensive than the numerical factorisation. It only
% depends on the sparsity pattern, but the pattern of the matrices solved in Ua changes from one Newton-Raphson iteration to the next:
% entries that are exactly zero are not stored in MATLAB sparse matrices, and the rows of constrained nodes change with the active set.
% Therefore a fixed pattern P is kept, which contains the pattern of all matrices of the same size and type solved so far. The values of
% A are copied into this pattern (by pardisoMex, with explicit zeros where A has no entry), and the analysis is only repeated if P has to
% be extended. P initially contains the pattern of A, the diagonal, and, if given, the pattern of Apattern (eg the full Jacobian, of
% which the rows of constrained nodes are removed in A). A separate PARDISO state, with its own pattern, is kept for each matrix size and
% type (eg uv and uvh systems); at most 6 states are kept, and the least recently used one is released when a new one is needed.
%
% (17 Oct 2026) Structural pattern: if StructuralPattern is given, it is used as the fixed pattern instead of the learned one, so that
% normally only one analysis is needed per mesh and type of system. StructuralPattern is either an assembly pattern of the mesh, as
% stored in MUA.uvAssemblyPattern or MUA.uvhAssemblyPattern (fields i0, j0 and neq, giving all entries the assembled matrices can have),
% or a structure with fields Pattern (such an assembly pattern) and Index (the unknowns of a reduced system), in which case the pattern
% of the reduced system is Pattern(Index,Index). It is ignored if it does not match the size of A, or if
% CtrlVar.LinearSolver.PARDISO.StructuralPattern=false. If A has an entry outside the structural pattern, the pattern is extended as
% described above, so that a wrong or outdated structural pattern costs time, but does not affect the solution.
%
% The relative residual of the solution is checked. If it is larger than CtrlVar.LinearSolver.PARDISO.ResidualTolerance, the analysis
% is repeated with the current values of A (the scaling and matching used for pivoting are computed in the analysis, from the values of
% A at that time), and the system solved again. If the residual is still too large, the system is solved with decomposition instead.
%
% The MKL runtime library (mkl_rt) is not part of MATLAB, and must be installed separately, eg from PyPI (pip install mkl). Its path is
% given by CtrlVar.LinearSolver.PARDISO.Library or, if empty, by the environment variable UA_MKL_RT. Within MATLAB, the TBB threading
% layer of MKL must be used (default), or the sequential layer; the OpenMP layer can not be used (see pardisoMex.c).
%
%%

persistent States Stats Clock

MaxStates=6 ;

if ischar(CtrlVar) || isstring(CtrlVar)
    switch string(CtrlVar)
        case "stats"
            if isempty(Stats) ; Stats=NewStats() ; end
            x=Stats ;
            x.States=numel(States) ;
        case "reset"
            for k=1:numel(States)
                try
                    pardisoMex('free',States(k).h) ;
                catch
                end
            end
            States=[] ; Stats=NewStats() ; x=[] ; StructuralPatternOf("reset") ;
        otherwise
            error("UaPardisoSolve:Usage","Unknown command %s.",CtrlVar)
    end
    return
end

if isempty(Stats) ; Stats=NewStats() ; end
if isempty(Clock) ; Clock=0 ; end
if nargin<4 ; Apattern=[] ; end
if nargin<5 || isempty(MatrixType) ; MatrixType="nonsymmetric" ; end
if nargin<6 ; StructuralPattern=[] ; end
switch lower(string(MatrixType))
    case "nonsymmetric" ; mtype=11 ;
    case "spd"          ; mtype=2 ;
    case "symmetric"    ; mtype=-2 ;
    otherwise
        error("UaPardisoSolve:MatrixType","MatrixType must be ""nonsymmetric"", ""spd"" or ""symmetric"".")
end
Opt=CtrlVar.LinearSolver.PARDISO ;
MinimumSize=20000 ; if isfield(Opt,"MinimumSize") && ~isempty(Opt.MinimumSize) ; MinimumSize=Opt.MinimumSize ; end
if issparse(b) ; b=full(b) ; end
n=size(A,1) ;

% small and diagonal systems
if isdiag(A)
    x=b./full(diag(A)) ;
    Stats.DiagonalSystems=Stats.DiagonalSystems+1 ;
    return
end
if n<MinimumSize
    x=A\b ;
    Stats.SmallSystems=Stats.SmallSystems+1 ;
    return
end

Lib=string(Opt.Library) ;
if Lib=="" ; Lib=string(getenv("UA_MKL_RT")) ; end
if Lib=="" || ~isfile(Lib)
    error("UaPardisoSolve:Library",...
        "The MKL runtime library (mkl_rt) was not found. Set CtrlVar.LinearSolver.PARDISO.Library, or the environment variable UA_MKL_RT, to its full path (eg ...\\Library\\bin\\mkl_rt.3.dll on Windows, .../lib/libmkl_rt.so.3 on Linux).")
end
pardisoMex('load',char(Lib),char(Opt.ThreadingLayer)) ;   % returns immediately if already loaded

Symmetric = mtype~=11 ;
if Symmetric
    Af=tril(A) ;          % the upper triangle in the CSR view of PARDISO
else
    Af=A ;
end

% (17 Oct 2026) the structural pattern of the mesh, if given
Pst=[] ; Sid=0 ;
UseStructural=~isfield(Opt,"StructuralPattern") || isempty(Opt.StructuralPattern) || Opt.StructuralPattern ;
if UseStructural && ~isempty(StructuralPattern)
    [Pst,Sid]=StructuralPatternOf(StructuralPattern,n,Symmetric) ;   % Sid=0 if it does not match
end

% the PARDISO state and the fixed pattern for matrices of this size and type
k=[] ;
if ~isempty(States) ; k=find([States.n]==n & [States.mtype]==mtype,1) ; end
if isempty(k)
    if numel(States)>=MaxStates
        [~,iLRU]=min([States.LastUsed]) ;
        try
            pardisoMex('free',States(iLRU).h) ;
        catch
        end
        States(iLRU)=[] ;
        Stats.StatesReleased=Stats.StatesReleased+1 ;
    end
    if Sid>0
        P=Pst | (Af~=0) ;
        Stats.StructuralPatterns=Stats.StructuralPatterns+1 ;
    else
        P=(Af~=0) | logical(speye(n)) ;
        if ~isempty(Apattern) ; P=P | PatternOf(Apattern,Symmetric) ; end
    end
    New=struct("n",n,"mtype",mtype,"h",pardisoMex('init',mtype),"P",P,"LastUsed",0,"Sid",Sid) ;
    if isempty(States) ; States=New ; else ; States(end+1)=New ; end
    k=numel(States) ;
elseif Sid>0 && States(k).Sid~=Sid
    % a (new) structural pattern for this size and type, eg after remeshing: it replaces the current pattern
    States(k).P=Pst | (Af~=0) ;
    States(k).Sid=Sid ;
    Stats.StructuralPatterns=Stats.StructuralPatterns+1 ;
end
Clock=Clock+1 ; States(k).LastUsed=Clock ;

% factorise, extending the pattern if A has entries outside it
try
    try
        r=pardisoMex('factorp',States(k).h,Af,States(k).P) ;
    catch ME
        if ME.identifier~="pardisoMex:handle" ; rethrow(ME) ; end
        % the solver states of pardisoMex have been released (eg pardisoMex was cleared from memory): start again
        States=struct("n",n,"mtype",mtype,"h",pardisoMex('init',mtype),"P",States(k).P,"LastUsed",Clock,"Sid",States(k).Sid) ;
        k=1 ;
        r=pardisoMex('factorp',States(k).h,Af,States(k).P) ;
    end
    if r<0
        States(k).P=States(k).P | (Af~=0) ;
        if ~isempty(Apattern) ; States(k).P=States(k).P | PatternOf(Apattern,Symmetric) ; end
        Stats.PatternExtensions=Stats.PatternExtensions+1 ;
        r=pardisoMex('factorp',States(k).h,Af,States(k).P) ;
    end
catch ME
    if mtype==2 && ME.identifier=="pardisoMex:pardiso"
        % Cholesky factorisation failed, presumably as the matrix is not positive definite: solve as symmetric indefinite
        warning("UaPardisoSolve:NotPositiveDefinite","PARDISO Cholesky factorisation failed (%s); solving the system as symmetric indefinite.",ME.message)
        Stats.SPDFallbacks=Stats.SPDFallbacks+1 ;
        x=UaPardisoSolve(CtrlVar,A,b,Apattern,"symmetric",StructuralPattern) ;
        return
    end
    rethrow(ME)
end
Stats.Analyses=Stats.Analyses+r ;
Stats.Factorisations=Stats.Factorisations+1 ;
x=pardisoMex('solvep',States(k).h,b) ;
Stats.Solves=Stats.Solves+1 ;

nb=max(norm(b,"fro"),realmin) ;
res=norm(A*x-b,"fro")/nb ;
if ~(res<=Opt.ResidualTolerance)
    % repeat the analysis with the current values of A, and solve again
    pardisoMex('factorp',States(k).h,Af,States(k).P,'reanalyse') ;
    x=pardisoMex('solvep',States(k).h,b) ;
    Stats.Analyses=Stats.Analyses+1 ;
    Stats.Factorisations=Stats.Factorisations+1 ;
    Stats.Resolves=Stats.Resolves+1 ;
    res=norm(A*x-b,"fro")/nb ;
    if ~(res<=Opt.ResidualTolerance)
        warning("UaPardisoSolve:Residual","PARDISO relative residual %g is above the tolerance %g; using decomposition instead.",...
            res,Opt.ResidualTolerance)
        x=decomposition(A)\b ;
        Stats.Fallbacks=Stats.Fallbacks+1 ;
        res=norm(A*x-b,"fro")/nb ;
    end
end
Stats.MaxResidual=max(Stats.MaxResidual,res) ;

end

%%
function P=PatternOf(M,Symmetric)
P=(M~=0) ;
if Symmetric ; P=tril(P) ; end
end

%%
function S=NewStats()
S=struct("Analyses",0,"Factorisations",0,"Solves",0,"PatternExtensions",0,"Resolves",0,"Fallbacks",0,"SPDFallbacks",0,...
    "SmallSystems",0,"DiagonalSystems",0,"StatesReleased",0,"StructuralPatterns",0,"MaxResidual",0,"States",0) ;
end

%%
function [P,id]=StructuralPatternOf(SP,n,Symmetric)

% (17 Oct 2026) Logical sparsity pattern (tril for the symmetric types, with the diagonal) from an assembly pattern of the mesh, possibly
% for a reduced system (SP.Pattern, SP.Index). Returns P=[] and id=0 if SP does not match the size n. The patterns are cached, and id
% identifies a pattern (the same id for the same structure, index set and type), so that a state only needs a new analysis when it changes.

persistent Cache NextId
if isempty(NextId) ; NextId=0 ; end
P=[] ; id=0 ;

if (ischar(SP) || isstring(SP)) && string(SP)=="reset"
    Cache=[] ;
    return
end

if isstruct(SP) && isfield(SP,"Pattern")
    Base=SP.Pattern ; Index=SP.Index(:) ;
else
    Base=SP ; Index=[] ;
end
if ~isstruct(Base) || ~all(isfield(Base,["i0","j0","neq"])) || isempty(Base.i0)
    return
end
if isempty(Index) ; m=Base.neq ; else ; m=numel(Index) ; end
if m~=n
    return
end

for c=1:numel(Cache)
    E=Cache(c) ;
    if E.neq==Base.neq && E.nk==numel(Base.i0) && E.Symmetric==Symmetric && isequal(E.Index,Index) ...
            && isequal(E.i0,Base.i0) && isequal(E.j0,Base.j0)
        P=E.P ; id=E.id ;
        return
    end
end

S=sparse(double(Base.i0),double(Base.j0),true,Base.neq,Base.neq) ;
if ~isempty(Index) ; S=S(Index,Index) ; end
if Symmetric ; S=tril(S) ; end
S=S | logical(speye(n)) ;

NextId=NextId+1 ;
E=struct("neq",Base.neq,"nk",numel(Base.i0),"Symmetric",Symmetric,"Index",Index,"i0",Base.i0,"j0",Base.j0,"P",S,"id",NextId) ;
if isempty(Cache) ; Cache=E ; else ; Cache=[Cache(max(1,end-2):end) E] ; end   % at most four patterns
P=S ; id=NextId ;

end
