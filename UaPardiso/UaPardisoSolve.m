function x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType,StructuralPattern,FactorIdToKeep)

%% Solves A x = b with the Intel oneMKL PARDISO sparse direct solver (through pardisoMex)
%
%   x=UaPardisoSolve(CtrlVar,A,b)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType)
%   x=UaPardisoSolve(CtrlVar,A,b,Apattern,MatrixType,StructuralPattern)
%   Stats=UaPardisoSolve("stats")     % numbers of analyses, factorisations, solves, pattern extensions, repeated solves, fallbacks, ...
%   UaPardisoSolve("reset")           % release the memory used by PARDISO and reset the statistics
%
%   Fac=UaPardisoSolve("factor",CtrlVar,A,Apattern,MatrixType,StructuralPattern)   % factorise A only (16 Nov 2026)
%   x=UaPardisoSolve("solve",Fac,b)                                                % solve A x = b with the factorisation Fac
%
% (16 Nov 2026) Reusing the factorisation: "factor" does the analysis (if needed) and the numerical factorisation, and returns a small
% structure Fac identifying the factorisation; "solve" then only carries out the solve phase, which is much cheaper than the
% factorisation. This is used when several systems with the same matrix are solved, eg the linearised forward and adjoint solves of
% the Hessian-vector products (HessianVectorProduct.m) and the modified Newton-Raphson method (reuse of dA in solveKApe.m). As one
% PARDISO state is kept for each matrix size and type, a factorisation can be replaced by that of another matrix of the same size and
% type in the meantime (eg a forward uv solve). Each factorisation is therefore given a unique number, and "solve" checks that the
% factorisation identified by Fac is still the current one; if not, A is factorised again (counted in the statistics as
% FactorReuseFailures), so that a replaced factorisation is never used. (24 Nov 2026) The new factorisation keeps the number of Fac, so
% that the following solves with Fac reuse it (previously every following solve with Fac factorised A again). The backward error of the
% solution is checked as for the other solves (see below).
% For small (n<MinimumSize) and diagonal systems Fac contains a MATLAB decomposition or the diagonal.
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
% (24 Nov 2026) Accuracy check: the normwise backward error of the solution,
%
%   eta = max_j ||A x_j - b_j||_inf / ( ||A||_inf ||x_j||_inf + ||b_j||_inf ) ,
%
% is compared with CtrlVar.LinearSolver.PARDISO.BackwardErrorTolerance. If eta is larger, and the analysis depends on the values of A
% (scaling and matching, iparm(11) and iparm(13); by default only for the nonsymmetric type), the analysis is repeated with the current
% values of A and the system solved again. For the symmetric types the analysis only depends on the pattern, so that repeating it would
% give the same factorisation. If eta is still too large, the system is solved with decomposition instead. For a factorisation used
% with "solve", this decomposition is kept (at most two), and the following solves with that factorisation use it directly.
%
% Previously the relative residual ||A x-b||/||b|| was checked (CtrlVar.LinearSolver.PARDISO.ResidualTolerance=1e-10, no longer used).
% That test is not scale invariant: ||A x-b||/||b|| can be as large as eta ||A|| ||x||/||b||, so that for badly scaled systems no
% backward stable solver can pass it. Eg in a logA, logC inversion (Amundsen Sea Embayment, 131k nodes, m=3, n=4) the SPD matrix K of
% the Hessian-vector products had diagonal values from 1.5e3 to 5e12 and ||A|| ||x||/||b||=4e9. PARDISO gave relative residuals of 1e-10
% to 1e-9, with backward errors of about 2e-19, and decomposition, used as fallback, a larger relative residual (2.6e-10) than PARDISO,
% so that the fallbacks (most Hessian-vector products) cost time without making the solutions more accurate.
%
% The MKL runtime library (mkl_rt) is not part of MATLAB, and must be installed separately, eg from PyPI (pip install mkl). Its path is
% given by CtrlVar.LinearSolver.PARDISO.Library or, if empty, by the environment variable UA_MKL_RT. Within MATLAB, the TBB threading
% layer of MKL must be used (default), or the sequential layer; the OpenMP layer can not be used (see pardisoMex.c).
%
%%

persistent States Stats Clock FactorCount FallbackCache

MaxStates=6 ;

Mode="factorise-and-solve" ;
if isempty(FactorCount) ; FactorCount=0 ; end

if ischar(CtrlVar) || isstring(CtrlVar)
    switch string(CtrlVar)
        case "factor"
            % Fac=UaPardisoSolve("factor",CtrlVar,A,Apattern,MatrixType,StructuralPattern): shift the inputs, and continue below
            Mode="factor" ;
            CtrlVar=A ; A=b ; b=zeros(size(A,1),0) ;
        case "solve"
            % x=UaPardisoSolve("solve",Fac,b)
            Fac=A ;
            if isempty(Stats) ; Stats=NewStats() ; end
            if isempty(Clock) ; Clock=0 ; end
            if issparse(b) ; b=full(b) ; end
            switch Fac.Type
                case "diagonal"
                    x=b./Fac.d ;
                case "decomposition"
                    x=Fac.dA\b ;
                case "PARDISO"
                    % (24 Nov 2026) a decomposition used instead of this factorisation, as its backward error was too large
                    c=[] ;
                    if ~isempty(FallbackCache) ; c=find([FallbackCache.FactorId]==Fac.FactorId,1) ; end
                    if ~isempty(c)
                        x=FallbackCache(c).dA\b ;
                        Stats.FallbackReuses=Stats.FallbackReuses+1 ;
                        return
                    end
                    k=[] ;
                    if ~isempty(States) ; k=find([States.FactorId]==Fac.FactorId,1) ; end
                    if isempty(k)
                        % the factorisation has been replaced in the meantime: factorise again (and solve). (24 Nov 2026) The new
                        % factorisation is of the same matrix, and keeps the number of Fac, so that the following solves with Fac reuse it.
                        Stats.FactorReuseFailures=Stats.FactorReuseFailures+1 ;
                        x=UaPardisoSolve(Fac.CtrlVar,Fac.A,b,Fac.Apattern,Fac.MatrixType,Fac.StructuralPattern,Fac.FactorId) ;
                        return
                    end
                    Clock=Clock+1 ; States(k).LastUsed=Clock ;
                    x=pardisoMex('solvep',States(k).h,b) ;
                    Stats.Solves=Stats.Solves+1 ;
                    Stats.FactorReuses=Stats.FactorReuses+1 ;
                    % (24 Nov 2026) backward error instead of relative residual; no new factorisation of the same matrix (see above)
                    Tol=BackwardErrorTolerance(Fac.CtrlVar.LinearSolver.PARDISO) ;
                    if isfield(Fac,"NormA") ; NormA=Fac.NormA ; else ; NormA=norm(Fac.A,inf) ; end
                    eta=BackwardError(Fac.A,x,b,NormA) ;
                    if ~(eta<=Tol) && AnalysisDependsOnValues(States(k).h)
                        % repeat the analysis with the values of this matrix, and solve again (still a factorisation of Fac.A: its number is kept)
                        if States(k).mtype==11 ; AfFac=Fac.A ; else ; AfFac=tril(Fac.A) ; end
                        pardisoMex('factorp',States(k).h,AfFac,States(k).P,'reanalyse') ;
                        x=pardisoMex('solvep',States(k).h,b) ;
                        Stats.Analyses=Stats.Analyses+1 ;
                        Stats.Factorisations=Stats.Factorisations+1 ;
                        Stats.Resolves=Stats.Resolves+1 ;
                        eta=BackwardError(Fac.A,x,b,NormA) ;
                    end
                    if ~(eta<=Tol)
                        fprintf("UaPardisoSolve: PARDISO backward error %g is above the tolerance %g; using decomposition instead.\n",eta,Tol)
                        dA=decomposition(Fac.A) ;
                        x=dA\b ;
                        FallbackCache=AddToFallbackCache(FallbackCache,Fac.FactorId,dA) ;
                        Stats.Fallbacks=Stats.Fallbacks+1 ;
                        eta=BackwardError(Fac.A,x,b,NormA) ;
                    end
                    Stats.MaxBackwardError=max(Stats.MaxBackwardError,eta) ;
                otherwise
                    error("UaPardisoSolve:Fac","Unknown factorisation type.")
            end
            return
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
            States=[] ; Stats=NewStats() ; x=[] ; StructuralPatternOf("reset") ; FallbackCache=[] ;
        otherwise
            error("UaPardisoSolve:Usage","Unknown command %s.",CtrlVar)
    end
    if Mode~="factor" ; return ; end
end

if isempty(Stats) ; Stats=NewStats() ; end
if isempty(Clock) ; Clock=0 ; end
if nargin<4 ; Apattern=[] ; end
if nargin<5 || isempty(MatrixType) ; MatrixType="nonsymmetric" ; end
if nargin<6 ; StructuralPattern=[] ; end
if nargin<7 || isempty(FactorIdToKeep) ; FactorIdToKeep=0 ; end   % (24 Nov 2026) see the "solve" mode
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
    if Mode=="factor" ; x=struct("Type","diagonal","d",full(diag(A))) ; return ; end
    x=b./full(diag(A)) ;
    Stats.DiagonalSystems=Stats.DiagonalSystems+1 ;
    return
end
if n<MinimumSize
    if Mode=="factor" ; x=struct("Type","decomposition","dA",decomposition(A)) ; return ; end
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
        P=Pst | StoredPattern(Af) ;
        Stats.StructuralPatterns=Stats.StructuralPatterns+1 ;
    else
        P=StoredPattern(Af) | logical(speye(n)) ;
        if ~isempty(Apattern) ; P=P | PatternOf(Apattern,Symmetric) ; end
    end
    New=struct("n",n,"mtype",mtype,"h",pardisoMex('init',mtype),"P",P,"LastUsed",0,"Sid",Sid,"FactorId",0) ;
    if isempty(States) ; States=New ; else ; States(end+1)=New ; end
    k=numel(States) ;
elseif Sid>0 && States(k).Sid~=Sid
    % a (new) structural pattern for this size and type, eg after remeshing: it replaces the current pattern
    States(k).P=Pst | StoredPattern(Af) ;
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
        States=struct("n",n,"mtype",mtype,"h",pardisoMex('init',mtype),"P",States(k).P,"LastUsed",Clock,"Sid",States(k).Sid,"FactorId",0) ;
        k=1 ;
        r=pardisoMex('factorp',States(k).h,Af,States(k).P) ;
    end
    if r<0
        States(k).P=States(k).P | StoredPattern(Af) ;
        if ~isempty(Apattern) ; States(k).P=States(k).P | PatternOf(Apattern,Symmetric) ; end
        Stats.PatternExtensions=Stats.PatternExtensions+1 ;
        r=pardisoMex('factorp',States(k).h,Af,States(k).P) ;
        if r<0
            error("UaPardisoSolve:Pattern","The matrix has entries outside the extended pattern (not factorised).")
        end
    end
catch ME
    if mtype==2 && ME.identifier=="pardisoMex:pardiso"
        % Cholesky factorisation failed, presumably as the matrix is not positive definite: solve as symmetric indefinite
        warning("UaPardisoSolve:NotPositiveDefinite","PARDISO Cholesky factorisation failed (%s); solving the system as symmetric indefinite.",ME.message)
        Stats.SPDFallbacks=Stats.SPDFallbacks+1 ;
        if Mode=="factor"   % (16 Nov 2026) return a factorisation (symmetric indefinite), not a solution
            x=UaPardisoSolve("factor",CtrlVar,A,Apattern,"symmetric",StructuralPattern) ;
        else
            x=UaPardisoSolve(CtrlVar,A,b,Apattern,"symmetric",StructuralPattern,FactorIdToKeep) ;   % (24 Nov 2026) FactorIdToKeep passed on
        end
        return
    end
    rethrow(ME)
end
Stats.Analyses=Stats.Analyses+r ;
Stats.Factorisations=Stats.Factorisations+1 ;
if FactorIdToKeep>0
    States(k).FactorId=FactorIdToKeep ;   % (24 Nov 2026) a new factorisation of the matrix of a given factorisation ("solve" mode)
else
    FactorCount=FactorCount+1 ;
    States(k).FactorId=FactorCount ;   % (16 Nov 2026) identifies this factorisation
end

if Mode=="factor"
    x=struct("Type","PARDISO","FactorId",States(k).FactorId,"A",A,"NormA",norm(A,inf),"CtrlVar",CtrlVar,"Apattern",Apattern,...
        "MatrixType",MatrixType,"StructuralPattern",StructuralPattern) ;
    return
end

x=pardisoMex('solvep',States(k).h,b) ;
Stats.Solves=Stats.Solves+1 ;

% (24 Nov 2026) Accuracy check with the normwise backward error (see above), instead of the relative residual
Tol=BackwardErrorTolerance(Opt) ;
NormA=norm(A,inf) ;
eta=BackwardError(A,x,b,NormA) ;
if ~(eta<=Tol) && AnalysisDependsOnValues(States(k).h)
    % Repeat the analysis with the current values of A (the scaling and matching used for pivoting are computed in the analysis, from
    % the values of A at that time), and solve again. (24 Nov 2026) Only if the analysis depends on the values of A, ie not for the
    % symmetric types with the default parameters, where it would give the same factorisation. The factorisation is of the same
    % matrix, so it keeps its number (previously a new number was given, which made a factorisation passed to "solve" mode stale).
    pardisoMex('factorp',States(k).h,Af,States(k).P,'reanalyse') ;
    x=pardisoMex('solvep',States(k).h,b) ;
    Stats.Analyses=Stats.Analyses+1 ;
    Stats.Factorisations=Stats.Factorisations+1 ;
    Stats.Resolves=Stats.Resolves+1 ;
    eta=BackwardError(A,x,b,NormA) ;
end
if ~(eta<=Tol)
    fprintf("UaPardisoSolve: PARDISO backward error %g is above the tolerance %g; using decomposition instead.\n",eta,Tol)
    dA=decomposition(A) ;
    x=dA\b ;
    if FactorIdToKeep>0
        % (24 Nov 2026) called from "solve" mode: the following solves with that factorisation use this decomposition
        FallbackCache=AddToFallbackCache(FallbackCache,FactorIdToKeep,dA) ;
    end
    Stats.Fallbacks=Stats.Fallbacks+1 ;
    eta=BackwardError(A,x,b,NormA) ;
end
Stats.MaxBackwardError=max(Stats.MaxBackwardError,eta) ;

end

%%
function P=PatternOf(M,Symmetric)
P=StoredPattern(M) ;
if Symmetric ; P=tril(P) ; end
end

%%
function S=NewStats()
S=struct("Analyses",0,"Factorisations",0,"Solves",0,"PatternExtensions",0,"Resolves",0,"Fallbacks",0,"SPDFallbacks",0,...
    "SmallSystems",0,"DiagonalSystems",0,"StatesReleased",0,"StructuralPatterns",0,"MaxBackwardError",0,"States",0,"FactorReuses",0,"FactorReuseFailures",0,...
    "FallbackReuses",0) ;   % (24 Nov 2026) MaxBackwardError replaces MaxResidual; FallbackReuses: solves with a kept decomposition
end

%%
function eta=BackwardError(A,x,b,NormA)

% (24 Nov 2026) Normwise backward error of the solution x of A x = b, the largest over the columns of b:
%   eta = max_j ||A x_j - b_j||_inf / ( ||A||_inf ||x_j||_inf + ||b_j||_inf )
% NormA=norm(A,inf). It is scale invariant, and of the order of the machine precision for a backward stable solver, whatever the
% conditioning and scaling of A. Non-finite values in x give eta=NaN, which fails the test.

if isempty(b)
    eta=0 ;
    return
end
r=A*x-b ;
eta=max( vecnorm(r,Inf,1) ./ max(NormA*vecnorm(x,Inf,1)+vecnorm(b,Inf,1),realmin) ) ;

end

%%
function Tol=BackwardErrorTolerance(Opt)

% (24 Nov 2026) CtrlVar.LinearSolver.PARDISO.BackwardErrorTolerance, with its default if not given. The former
% CtrlVar.LinearSolver.PARDISO.ResidualTolerance (a relative residual) is no longer used; a note is printed once if it is given.

persistent Noted

Tol=1e-12 ;
if isfield(Opt,"BackwardErrorTolerance") && ~isempty(Opt.BackwardErrorTolerance)
    Tol=Opt.BackwardErrorTolerance ;
end
if isempty(Noted) && isfield(Opt,"ResidualTolerance") && ~isempty(Opt.ResidualTolerance)
    Noted=true ;
    fprintf("UaPardisoSolve: CtrlVar.LinearSolver.PARDISO.ResidualTolerance is no longer used; the accuracy check now uses the backward error, with CtrlVar.LinearSolver.PARDISO.BackwardErrorTolerance=%g.\n",Tol)
end

end

%%
function tf=AnalysisDependsOnValues(h)

% (24 Nov 2026) True if the analysis (phase 11) of the PARDISO state h uses the values of the matrix, ie scaling (iparm(11)) or
% weighted matching (iparm(13)) is switched on. By default this is the case for the nonsymmetric type (11) only; for the symmetric
% types the analysis then only depends on the pattern, and repeating it gives the same factorisation.

ip=pardisoMex('iparm',h) ;
tf = ip(11)~=0 || ip(13)~=0 ;

end

%%
function C=AddToFallbackCache(C,FactorId,dA)

% (24 Nov 2026) Keeps the decomposition dA used instead of the PARDISO factorisation with the number FactorId (at most two are kept, as
% these can be large).

E=struct("FactorId",FactorId,"dA",dA) ;
if isempty(C)
    C=E ;
else
    C=[C(end) E] ;
end

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
    if E.neq==Base.neq && E.nk==numel(Base.i0) && E.Symmetric==Symmetric && isequal(E.Index,Index)
        if isfield(Base,"TopologyKey") && ~isempty(E.TopologyKey)
            Same=isequal(E.TopologyKey,Base.TopologyKey) ;   % (17 Oct 2026) key of the mesh (see MeshKey.m)
        else
            Same=isequal(E.i0,Base.i0) && isequal(E.j0,Base.j0) ;
        end
        if Same
            P=E.P ; id=E.id ;
            return
        end
    end
end

S=sparse(double(Base.i0),double(Base.j0),true,Base.neq,Base.neq) ;
if ~isempty(Index) ; S=S(Index,Index) ; end
if Symmetric ; S=tril(S) ; end
S=S | logical(speye(n)) ;

NextId=NextId+1 ;
TK=[] ; if isfield(Base,"TopologyKey") ; TK=Base.TopologyKey ; end
E=struct("neq",Base.neq,"nk",numel(Base.i0),"Symmetric",Symmetric,"Index",Index,"i0",Base.i0,"j0",Base.j0,"P",S,"id",NextId,"TopologyKey",TK) ;
if isempty(Cache) ; Cache=E ; else ; Cache=[Cache(max(1,end-2):end) E] ; end   % at most four patterns
P=S ; id=NextId ;

end

%%
function P=StoredPattern(M)

% (16 Nov 2026) Pattern of the stored entries of the sparse matrix M, including entries that are stored with the value zero. pardisoMex
% copies all stored entries of A into the fixed pattern, so the pattern must contain the explicitly stored zeros as well (eg the
% assembled uv Jacobian can contain a few of these). Previously (M~=0) was used, which excludes them, so that for such a matrix the
% pattern could not be extended to contain A, and A was not factorised.

[i,j]=find(M) ;        % find also returns the stored zeros
P=sparse(i,j,true,size(M,1),size(M,2)) ;

end
