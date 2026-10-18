/*
 * pardisoMex.c  --  MATLAB interface to the Intel oneMKL PARDISO sparse direct solver
 *
 * (16 Oct 2026) Test version. The MKL runtime library (mkl_rt) is loaded at run time, so no import libraries are needed.
 * MATLAB stores sparse matrices in compressed sparse column (CSC) format with 64-bit indices. The CSC arrays of A are the
 * compressed sparse row (CSR) arrays of A', which are passed to the 64-bit interface of PARDISO (pardiso_64) without copying,
 * with 0-based indexing (iparm(35)=1), and PARDISO is asked to solve the transposed system (iparm(12)=2), ie A x = b.
 *
 * Usage:
 *   pardisoMex('load', mklrtPath)          % load mkl_rt (eg ...\Library\bin\mkl_rt.3.dll), with the TBB threading layer
 *   pardisoMex('load', mklrtPath, layer)   % layer: "TBB" (default, multithreaded) or "SEQUENTIAL" (one thread)
 *   layer = pardisoMex('layer')            % threading layer in use
 *   v = pardisoMex('version')              % MKL version string
 *   pardisoMex('setthreads', n)            % number of threads used by MKL
 *   n = pardisoMex('getthreads')
 *   h = pardisoMex('init', mtype)          % new solver state (mtype 11: real nonsymmetric), returns handle
 *   ip = pardisoMex('iparm', h)            % current iparm (1-based as in the MKL documentation)
 *   pardisoMex('setiparm', h, i, value)    % set iparm(i), 1-based
 *   pardisoMex('analyse', h, A)            % phase 11: reordering and symbolic factorisation
 *   pardisoMex('factor', h, A)             % phase 22: numerical factorisation (A must have the pattern used in the analysis)
 *   r = pardisoMex('factor', h, A, 'auto') % as above, but if the pattern of A differs from that of the analysis, the analysis is
 *                                          % repeated first; r=1 if this was done, otherwise r=0
 *   x = pardisoMex('solve', h, A, b)       % phase 33: solve A x = b (b may have several columns)
 *
 *   Fixed-pattern mode (a solver state is used either in this mode or in the one above, not both):
 *   r = pardisoMex('factorp', h, A, P)     % P: sparse matrix whose pattern contains that of A (values of P are not used). The
 *                                          % pattern of P is copied, the values of A are copied into it (with explicit zeros where A
 *                                          % has no entry), the analysis is done if P has changed since the last analysis, and A is
 *                                          % factorised. r=1: analysis done, r=0: analysis reused, r=-1: A has an entry outside the
 *                                          % pattern of P (nothing factorised; extend P and call again).
 *   r = pardisoMex('factorp', h, A, P, 'reanalyse')   % as above, but the analysis is always repeated (with the current values)
 *   x = pardisoMex('solvep', h, b)         % solve with the matrix of the last 'factorp'
 *   pardisoMex('free', h)                  % phase -1: release all memory of this solver state
 *
 * Threading layer: MKL's OpenMP layer (mkl_intel_thread) must not be used within MATLAB: it requires a version of the OpenMP
 * runtime (libiomp5md.dll) that may differ from the one MATLAB has already loaded, and MKL then aborts the process. The TBB
 * layer uses the TBB library already loaded by MATLAB (tested with R2026b and MKL 2026.1). The layer is set by this interface,
 * before the first computation, with MKL_Set_Threading_Layer. With the TBB layer, MKL uses all cores; 'setthreads' has no effect.
 *
 * Compile (MinGW64):  mex -R2018a pardisoMex.c
 *
 * Linux (16 Oct 2026, NOT YET TESTED WITHIN MATLAB; the loading and the PARDISO calls were tested outside MATLAB with MKL 2026.2):
 *   The library is libmkl_rt.so.N (N=2 for MKL 2021-2025, N=3 for MKL 2026), eg .../lib/libmkl_rt.so.3 after pip install mkl. It is
 *   loaded with dlopen(RTLD_NOW|RTLD_LOCAL|RTLD_DEEPBIND), so that MKL prefers its own symbols to those of the (reduced) MKL that MATLAB
 *   has loaded. libmkl_rt loads its other libraries at run time by name; as they have no RUNPATH, and LD_LIBRARY_PATH can not be changed
 *   within MATLAB, libmkl_core, libmkl_intel_lp64/ilp64 and the threading layer are loaded first, from the folder of libmkl_rt. The TBB
 *   layer needs libtbb.so.12, normally that of MATLAB; if that is not found, the one in the folder of MKL is used. As on Windows the
 *   OpenMP layer is not used. If MATLAB crashes or MKL can not be loaded, try the SEQUENTIAL layer first.
 *   Compile with gcc:  mex -R2018a pardisoMex.c -ldl   (-ldl is only needed for glibc older than 2.34)
 */

#if !defined(_WIN32) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE            /* for RTLD_DEEPBIND (glibc) */
#endif
#include "mex.h"
#include <string.h>
#include <stdio.h>

/* ---- loading of the MKL runtime library at run time: Windows and Linux ---- */
#ifdef _WIN32

#include <windows.h>
typedef HMODULE LibHandle;
#define STRICMP _stricmp

static LibHandle libOpen(const char *path, int layer, char *msg, size_t len)
{
    LibHandle h;
    (void) layer;
    h = LoadLibraryExA(path, NULL, LOAD_WITH_ALTERED_SEARCH_PATH);   /* dependencies are searched in the folder of path */
    if (h == NULL) snprintf(msg, len, "Windows error %lu", GetLastError());
    return h;
}
static void *libSym(LibHandle h, const char *name) { return (void *) GetProcAddress(h, name); }
static void  libClose(LibHandle h) { FreeLibrary(h); }

#else   /* Linux (not tested within MATLAB, see the notes at the top) */

#include <dlfcn.h>
#include <strings.h>
typedef void *LibHandle;
#define STRICMP strcasecmp

static int dlFlags(void)
{
    int flags = RTLD_NOW | RTLD_LOCAL;
#ifdef RTLD_DEEPBIND
    flags |= RTLD_DEEPBIND;    /* MKL prefers its own symbols to those of the MKL already loaded by MATLAB */
#endif
    return flags;
}

static LibHandle libOpen(const char *path, int layer, char *msg, size_t len)
{
    /* path is eg .../lib/libmkl_rt.so.3. libmkl_rt loads the other MKL libraries at run time by name, and as these have no RUNPATH and
       LD_LIBRARY_PATH can not be changed within MATLAB, they are loaded here first, by full path from the same folder and with the same
       version suffix; libmkl_rt then finds them already loaded. Failures are ignored here (libmkl_rt then searches itself). */
    char dir[2048], suffix[64], dep[2300];
    const char *slash = strrchr(path, '/'), *base = slash ? slash + 1 : path, *p = strstr(base, ".so");
    LibHandle h;
    if (slash != NULL && p != NULL && (size_t) (slash - path) < sizeof(dir) && strlen(p) < sizeof(suffix)) {
        const char *deps[5];
        int k, nd = 0;
        memcpy(dir, path, (size_t) (slash - path)); dir[slash - path] = 0;
        strcpy(suffix, p);                                   /* eg ".so.3" */
        deps[nd++] = "libmkl_core";
        deps[nd++] = "libmkl_intel_lp64";
        deps[nd++] = "libmkl_intel_ilp64";
        deps[nd++] = (layer == 4) ? "libmkl_tbb_thread" : "libmkl_sequential";
        for (k = 0; k < nd; k++) {
            snprintf(dep, sizeof(dep), "%s/%s%s", dir, deps[k], suffix);
            if (dlopen(dep, dlFlags()) == NULL && layer == 4 && k == nd - 1) {
                /* libmkl_tbb_thread needs libtbb.so.12: normally that of MATLAB (already loaded, or in its library folder);
                   only if it is not found, the one in the folder of MKL is used */
                char tbb[2300];
                snprintf(tbb, sizeof(tbb), "%s/libtbb.so.12", dir);
                if (dlopen(tbb, RTLD_NOW | RTLD_GLOBAL) != NULL) dlopen(dep, dlFlags());
            }
        }
    }
    h = dlopen(path, dlFlags());
    if (h == NULL) snprintf(msg, len, "%s", dlerror());
    return h;
}
static void *libSym(LibHandle h, const char *name) { return dlsym(h, name); }
static void  libClose(LibHandle h) { dlclose(h); }

#endif

typedef long long int i64;

typedef void (*pardiso64_t)(void *pt, const i64 *maxfct, const i64 *mnum, const i64 *mtype, const i64 *phase, const i64 *n,
                            const void *a, const i64 *ia, const i64 *ja, i64 *perm, const i64 *nrhs, i64 *iparm,
                            const i64 *msglvl, void *b, void *x, i64 *error);
typedef void (*pardisoinit_t)(void *pt, const int *mtype, int *iparm);
typedef void (*setthreads_t)(int n);
typedef int  (*getthreads_t)(void);
typedef void (*version_t)(char *buf, int len);
typedef int  (*setlayer_t)(int layer);

#define MAXSOLVERS 16

typedef struct {
    int    used;
    void  *pt[64];
    i64    iparm[64];
    i64    mtype;
    i64    n;
    i64    nnz;
    unsigned long long hash;   /* checksum of the pattern (column pointers and row indices) used in the analysis */
    int    analysed;
    int    factorised;
    int    patternMode;        /* 1: fixed-pattern mode ('factorp', 'solvep') */
    i64   *pia, *pja;          /* fixed-pattern mode: copies of the column pointers and row indices of P (persistent memory) */
    double *pval;              /* fixed-pattern mode: values of A in the pattern of P (persistent memory) */
} SolverState;

static LibHandle     hMKL = NULL;
static pardiso64_t   f_pardiso64 = NULL;
static pardisoinit_t f_pardisoinit = NULL;
static setthreads_t  f_setthreads = NULL;
static getthreads_t  f_getthreads = NULL;
static version_t     f_version = NULL;
static setlayer_t    f_setlayer = NULL;
static int           threadingLayer = -1;   /* 1: sequential, 4: TBB (MKL_THREADING_SEQUENTIAL, MKL_THREADING_TBB) */
static SolverState   S[MAXSOLVERS];
static int           atExitRegistered = 0;

static unsigned long long patternHash(const mxArray *A)
{
    /* FNV-1a over the column pointers and row indices */
    const mwIndex *jc = mxGetJc(A), *ir = mxGetIr(A);
    mwSize n = mxGetN(A), k, nnz = jc[n];
    unsigned long long h = 1469598103934665603ULL;
    for (k = 0; k <= n; k++)  { h ^= (unsigned long long) jc[k]; h *= 1099511628211ULL; }
    for (k = 0; k < nnz; k++) { h ^= (unsigned long long) ir[k]; h *= 1099511628211ULL; }
    return h;
}

static void callPardisoRaw(SolverState *s, i64 phase, const double *aval, const i64 *ia, const i64 *ja, double *b, double *x, i64 nrhs)
{
    i64 maxfct = 1, mnum = 1, msglvl = 0, error = 0;
    i64 n = s->n, idum = 0;
    double ddum = 0.0;
    f_pardiso64(s->pt, &maxfct, &mnum, &s->mtype, &phase, &n, aval ? (const void *) aval : (const void *) &ddum,
                ia ? ia : &idum, ja ? ja : &idum, NULL, &nrhs, s->iparm, &msglvl,
                b ? (void *) b : (void *) &ddum, x ? (void *) x : (void *) &ddum, &error);
    if (error != 0) {
        mexErrMsgIdAndTxt("pardisoMex:pardiso", "PARDISO returned error %lld in phase %lld (see the MKL documentation).",
                          error, phase);
    }
}

static void freePatternBuffers(SolverState *s)
{
    if (s->pia)  { mxFree(s->pia);  s->pia = NULL; }
    if (s->pja)  { mxFree(s->pja);  s->pja = NULL; }
    if (s->pval) { mxFree(s->pval); s->pval = NULL; }
}

static void callPardiso(SolverState *s, i64 phase, const mxArray *A, double *b, double *x, i64 nrhs)
{
    i64 maxfct = 1, mnum = 1, msglvl = 0, error = 0;
    i64 n = s->n, idum = 0;
    double ddum = 0.0;
    const void *a = (const void *) &ddum;          /* valid dummy pointers for phases without a matrix (eg -1) */
    const i64 *ia = &idum, *ja = &idum;
    if (A != NULL) {
        a  = (const void *) mxGetDoubles(A);
        ia = (const i64 *) mxGetJc(A);
        ja = (const i64 *) mxGetIr(A);
    }
    f_pardiso64(s->pt, &maxfct, &mnum, &s->mtype, &phase, &n, a, ia, ja, NULL, &nrhs, s->iparm, &msglvl,
                b ? (void *) b : (void *) &ddum, x ? (void *) x : (void *) &ddum, &error);
    if (error != 0) {
        mexErrMsgIdAndTxt("pardisoMex:pardiso", "PARDISO returned error %lld in phase %lld (see the MKL documentation).",
                          error, phase);
    }
}

static void releaseAll(void)
{
    int k;
    for (k = 0; k < MAXSOLVERS; k++) {
        if (S[k].used && f_pardiso64 != NULL) {
            i64 maxfct = 1, mnum = 1, msglvl = 0, error = 0, phase = -1, nrhs = 1, n = S[k].n, idum = 0;
            double ddum = 0.0;
            if (!S[k].analysed) { freePatternBuffers(&S[k]); memset(&S[k], 0, sizeof(SolverState)); continue; }
            f_pardiso64(S[k].pt, &maxfct, &mnum, &S[k].mtype, &phase, &n, &ddum, &idum, &idum, NULL, &nrhs, S[k].iparm,
                        &msglvl, &ddum, &ddum, &error);
        }
        freePatternBuffers(&S[k]);
        memset(&S[k], 0, sizeof(SolverState));
    }
}

static SolverState *getSolver(const mxArray *h)
{
    int k;
    if (!mxIsNumeric(h) || mxGetNumberOfElements(h) != 1) mexErrMsgIdAndTxt("pardisoMex:handle", "Invalid handle.");
    k = (int) mxGetScalar(h) - 1;
    if (k < 0 || k >= MAXSOLVERS || !S[k].used) mexErrMsgIdAndTxt("pardisoMex:handle", "Invalid or released handle.");
    return &S[k];
}

static void checkMatrix(SolverState *s, const mxArray *A, int requireSamePattern)
{
    if (!mxIsSparse(A) || !mxIsDouble(A) || mxIsComplex(A))
        mexErrMsgIdAndTxt("pardisoMex:matrix", "A must be a real, double, sparse matrix.");
    if (mxGetM(A) != mxGetN(A)) mexErrMsgIdAndTxt("pardisoMex:matrix", "A must be square.");
    if (requireSamePattern) {
        i64 n = (i64) mxGetN(A), nnz = (i64) mxGetJc(A)[mxGetN(A)];
        if (!s->analysed) mexErrMsgIdAndTxt("pardisoMex:state", "Call 'analyse' first.");
        if (n != s->n || nnz != s->nnz)
            mexErrMsgIdAndTxt("pardisoMex:pattern", "A has a different size or number of nonzeros than in the analysis.");
    }
}

static void requireLoaded(void)
{
    if (hMKL == NULL || f_pardiso64 == NULL) mexErrMsgIdAndTxt("pardisoMex:notloaded", "Call pardisoMex('load', mklrtPath) first.");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char cmd[32];
    if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("pardisoMex:usage", "First argument must be a command string.");
    mxGetString(prhs[0], cmd, sizeof(cmd));
    if (!atExitRegistered) { mexAtExit(releaseAll); atExitRegistered = 1; }

    if (strcmp(cmd, "load") == 0) {
        char path[2048];
        char layerName[32];
        int layer = 4;                                   /* MKL_THREADING_TBB */
        if (nrhs < 2 || !mxIsChar(prhs[1])) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('load', mklrtPath [, layer])");
        mxGetString(prhs[1], path, sizeof(path));
        if (nrhs >= 3) {
            if (!mxIsChar(prhs[2])) mexErrMsgIdAndTxt("pardisoMex:usage", "layer must be \"TBB\" or \"SEQUENTIAL\".");
            mxGetString(prhs[2], layerName, sizeof(layerName));
            if (STRICMP(layerName, "TBB") == 0) layer = 4;
            else if (STRICMP(layerName, "SEQUENTIAL") == 0) layer = 1;
            else mexErrMsgIdAndTxt("pardisoMex:usage", "layer must be \"TBB\" or \"SEQUENTIAL\" (the OpenMP layer can not be used within MATLAB).");
        }
        if (hMKL != NULL && threadingLayer != layer)
            mexErrMsgIdAndTxt("pardisoMex:load", "MKL is already loaded with another threading layer; restart MATLAB to change it.");
        if (hMKL == NULL) {
            char msg[1024];
            msg[0] = 0;
            hMKL = libOpen(path, layer, msg, sizeof(msg));
            if (hMKL == NULL) mexErrMsgIdAndTxt("pardisoMex:load", "Could not load %s (%s).", path, msg);
            f_pardiso64   = (pardiso64_t)   libSym(hMKL, "pardiso_64");
            f_pardisoinit = (pardisoinit_t) libSym(hMKL, "pardisoinit");
            f_setthreads  = (setthreads_t)  libSym(hMKL, "MKL_Set_Num_Threads");
            f_getthreads  = (getthreads_t)  libSym(hMKL, "MKL_Get_Max_Threads");
            f_version     = (version_t)     libSym(hMKL, "MKL_Get_Version_String");
            f_setlayer    = (setlayer_t)    libSym(hMKL, "MKL_Set_Threading_Layer");
            if (f_pardiso64 == NULL || f_pardisoinit == NULL) {
                libClose(hMKL); hMKL = NULL; f_pardiso64 = NULL;
                mexErrMsgIdAndTxt("pardisoMex:load", "pardiso_64 or pardisoinit not found in %s.", path);
            }
            if (f_setlayer == NULL) {
                libClose(hMKL); hMKL = NULL; f_pardiso64 = NULL;
                mexErrMsgIdAndTxt("pardisoMex:load", "MKL_Set_Threading_Layer not found in %s.", path);
            }
            /* must be done before the first computation; returns the layer in use, or -1 */
            if (f_setlayer(layer) != layer) {
                mexErrMsgIdAndTxt("pardisoMex:load", "Could not set the MKL threading layer %d.", layer);
            }
            threadingLayer = layer;
        }
        return;
    }

    if (strcmp(cmd, "layer") == 0) {
        requireLoaded();
        plhs[0] = mxCreateString(threadingLayer == 4 ? "TBB" : (threadingLayer == 1 ? "SEQUENTIAL" : "unknown"));
        return;
    }

    if (strcmp(cmd, "version") == 0) {
        char buf[256];
        requireLoaded();
        buf[0] = 0;
        if (f_version) f_version(buf, (int) sizeof(buf));
        plhs[0] = mxCreateString(buf);
        return;
    }

    if (strcmp(cmd, "setthreads") == 0) {
        requireLoaded();
        if (nrhs < 2) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('setthreads', n)");
        if (f_setthreads) f_setthreads((int) mxGetScalar(prhs[1]));
        return;
    }

    if (strcmp(cmd, "getthreads") == 0) {
        requireLoaded();
        plhs[0] = mxCreateDoubleScalar(f_getthreads ? (double) f_getthreads() : -1.0);
        return;
    }

    if (strcmp(cmd, "init") == 0) {
        int k, mtype32, iparm32[64], i;
        requireLoaded();
        mtype32 = (nrhs >= 2) ? (int) mxGetScalar(prhs[1]) : 11;
        for (k = 0; k < MAXSOLVERS; k++) if (!S[k].used) break;
        if (k == MAXSOLVERS) mexErrMsgIdAndTxt("pardisoMex:init", "Too many solver states; release some with 'free'.");
        memset(&S[k], 0, sizeof(SolverState));
        memset(iparm32, 0, sizeof(iparm32));
        f_pardisoinit(S[k].pt, &mtype32, iparm32);   /* default parameters for this matrix type */
        for (i = 0; i < 64; i++) S[k].iparm[i] = (i64) iparm32[i];
        memset(S[k].pt, 0, sizeof(S[k].pt));         /* pt must be zero before the first call of pardiso_64 */
        S[k].iparm[0]  = 1;     /* non-default values supplied */
        S[k].iparm[11] = 2;     /* solve the transposed system: the CSC arrays of A are the CSR arrays of A' */
        S[k].iparm[34] = 1;     /* zero-based indexing */
        S[k].mtype = (i64) mtype32;
        S[k].used = 1;
        plhs[0] = mxCreateDoubleScalar((double) (k + 1));
        return;
    }

    if (strcmp(cmd, "iparm") == 0) {
        SolverState *s; int i; double *p;
        requireLoaded();
        if (nrhs < 2) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('iparm', h)");
        s = getSolver(prhs[1]);
        plhs[0] = mxCreateDoubleMatrix(64, 1, mxREAL);
        p = mxGetDoubles(plhs[0]);
        for (i = 0; i < 64; i++) p[i] = (double) s->iparm[i];
        return;
    }

    if (strcmp(cmd, "setiparm") == 0) {
        SolverState *s; int i;
        requireLoaded();
        if (nrhs < 4) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('setiparm', h, i, value)");
        s = getSolver(prhs[1]);
        i = (int) mxGetScalar(prhs[2]);
        if (i < 2 || i > 64) mexErrMsgIdAndTxt("pardisoMex:usage", "iparm index must be between 2 and 64 (1-based).");
        if (i == 12 || i == 35) mexErrMsgIdAndTxt("pardisoMex:usage", "iparm(12) and iparm(35) are fixed by this interface.");
        s->iparm[i - 1] = (i64) mxGetScalar(prhs[3]);
        return;
    }

    if (strcmp(cmd, "analyse") == 0) {
        SolverState *s;
        requireLoaded();
        if (nrhs < 3) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('analyse', h, A)");
        s = getSolver(prhs[1]);
        if (s->patternMode) mexErrMsgIdAndTxt("pardisoMex:mode", "This solver state is used in fixed-pattern mode ('factorp', 'solvep').");
        checkMatrix(s, prhs[2], 0);
        if (s->analysed) {   /* release the previous analysis and factorisation */
            callPardiso(s, -1, NULL, NULL, NULL, 1);
            memset(s->pt, 0, sizeof(s->pt));
            s->analysed = 0; s->factorised = 0;
        }
        s->n = (i64) mxGetN(prhs[2]);
        s->nnz = (i64) mxGetJc(prhs[2])[mxGetN(prhs[2])];
        s->hash = patternHash(prhs[2]);
        callPardiso(s, 11, prhs[2], NULL, NULL, 1);
        s->analysed = 1; s->factorised = 0;
        return;
    }

    if (strcmp(cmd, "factor") == 0) {
        SolverState *s;
        requireLoaded();
        int autoAnalyse = 0, reanalysed = 0;
        if (nrhs < 3) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('factor', h, A [, 'auto'])");
        s = getSolver(prhs[1]);
        if (s->patternMode) mexErrMsgIdAndTxt("pardisoMex:mode", "This solver state is used in fixed-pattern mode ('factorp', 'solvep').");
        if (nrhs >= 4) {
            char opt[16];
            if (!mxIsChar(prhs[3])) mexErrMsgIdAndTxt("pardisoMex:usage", "The fourth argument must be 'auto'.");
            mxGetString(prhs[3], opt, sizeof(opt));
            if (STRICMP(opt, "auto") != 0) mexErrMsgIdAndTxt("pardisoMex:usage", "The fourth argument must be 'auto'.");
            autoAnalyse = 1;
        }
        if (autoAnalyse) {
            unsigned long long h;
            checkMatrix(s, prhs[2], 0);
            h = patternHash(prhs[2]);
            if (!s->analysed || (i64) mxGetN(prhs[2]) != s->n || (i64) mxGetJc(prhs[2])[mxGetN(prhs[2])] != s->nnz || h != s->hash) {
                if (s->analysed) {
                    callPardiso(s, -1, NULL, NULL, NULL, 1);
                    memset(s->pt, 0, sizeof(s->pt));
                    s->analysed = 0; s->factorised = 0;
                }
                s->n = (i64) mxGetN(prhs[2]);
                s->nnz = (i64) mxGetJc(prhs[2])[mxGetN(prhs[2])];
                s->hash = h;
                callPardiso(s, 11, prhs[2], NULL, NULL, 1);
                s->analysed = 1;
                reanalysed = 1;
            }
        } else {
            checkMatrix(s, prhs[2], 1);
        }
        callPardiso(s, 22, prhs[2], NULL, NULL, 1);
        s->factorised = 1;
        if (nlhs > 0) plhs[0] = mxCreateDoubleScalar((double) reanalysed);
        return;
    }

    if (strcmp(cmd, "solve") == 0) {
        SolverState *s; i64 nr; double *bcopy, *x;
        requireLoaded();
        if (nrhs < 4) mexErrMsgIdAndTxt("pardisoMex:usage", "x = pardisoMex('solve', h, A, b)");
        s = getSolver(prhs[1]);
        if (s->patternMode) mexErrMsgIdAndTxt("pardisoMex:mode", "This solver state is used in fixed-pattern mode ('factorp', 'solvep').");
        checkMatrix(s, prhs[2], 1);
        if (!s->factorised) mexErrMsgIdAndTxt("pardisoMex:state", "Call 'factor' first.");
        if (mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (i64) mxGetM(prhs[3]) != s->n)
            mexErrMsgIdAndTxt("pardisoMex:rhs", "b must be a real, full, double matrix with as many rows as A.");
        nr = (i64) mxGetN(prhs[3]);
        bcopy = (double *) mxMalloc((size_t) (s->n * nr) * sizeof(double));
        memcpy(bcopy, mxGetDoubles(prhs[3]), (size_t) (s->n * nr) * sizeof(double));
        plhs[0] = mxCreateDoubleMatrix((mwSize) s->n, (mwSize) nr, mxREAL);
        x = mxGetDoubles(plhs[0]);
        callPardiso(s, 33, prhs[2], bcopy, x, nr);
        mxFree(bcopy);
        return;
    }

    if (strcmp(cmd, "free") == 0) {
        SolverState *s;
        requireLoaded();
        if (nrhs < 2) mexErrMsgIdAndTxt("pardisoMex:usage", "pardisoMex('free', h)");
        s = getSolver(prhs[1]);
        if (s->analysed) callPardiso(s, -1, NULL, NULL, NULL, 1);
        freePatternBuffers(s);
        memset(s, 0, sizeof(SolverState));
        return;
    }

    if (strcmp(cmd, "factorp") == 0) {
        SolverState *s;
        const mxArray *A, *P;
        const mwIndex *ajc, *air, *pjc, *pir;
        const double *aval;
        mwSize n, j;
        i64 pnnz;
        int force = 0, reanalysed = 0;
        unsigned long long h;
        requireLoaded();
        if (nrhs < 4) mexErrMsgIdAndTxt("pardisoMex:usage", "r = pardisoMex('factorp', h, A, P [, 'reanalyse'])");
        s = getSolver(prhs[1]);
        if (!s->patternMode && s->analysed) mexErrMsgIdAndTxt("pardisoMex:mode", "This solver state is used with 'analyse'/'factor'/'solve'.");
        A = prhs[2]; P = prhs[3];
        checkMatrix(s, A, 0);
        if (!mxIsSparse(P) || mxGetM(P) != mxGetM(A) || mxGetN(P) != mxGetN(A))
            mexErrMsgIdAndTxt("pardisoMex:pattern", "P must be a sparse matrix of the same size as A.");
        if (nrhs >= 5) {
            char opt[16];
            if (!mxIsChar(prhs[4])) mexErrMsgIdAndTxt("pardisoMex:usage", "The fifth argument must be 'reanalyse'.");
            mxGetString(prhs[4], opt, sizeof(opt));
            if (STRICMP(opt, "reanalyse") != 0) mexErrMsgIdAndTxt("pardisoMex:usage", "The fifth argument must be 'reanalyse'.");
            force = 1;
        }
        n = mxGetN(A);
        ajc = mxGetJc(A); air = mxGetIr(A); aval = mxGetDoubles(A);
        pjc = mxGetJc(P); pir = mxGetIr(P); pnnz = (i64) pjc[n];
        h = patternHash(P);

        /* new or changed pattern: copy it, and release the previous analysis */
        if (!s->patternMode || s->n != (i64) n || s->nnz != pnnz || s->hash != h) {
            if (s->analysed) { callPardisoRaw(s, -1, NULL, NULL, NULL, NULL, NULL, 1); memset(s->pt, 0, sizeof(s->pt)); }
            freePatternBuffers(s);
            s->analysed = 0; s->factorised = 0;
            s->pia  = (i64 *) mxMalloc((size_t) (n + 1) * sizeof(i64));
            s->pja  = (i64 *) mxMalloc((size_t) (pnnz > 0 ? pnnz : 1) * sizeof(i64));
            s->pval = (double *) mxMalloc((size_t) (pnnz > 0 ? pnnz : 1) * sizeof(double));
            mexMakeMemoryPersistent(s->pia); mexMakeMemoryPersistent(s->pja); mexMakeMemoryPersistent(s->pval);
            for (j = 0; j <= n; j++) s->pia[j] = (i64) pjc[j];
            for (j = 0; j < (mwSize) pnnz; j++) s->pja[j] = (i64) pir[j];
            s->n = (i64) n; s->nnz = pnnz; s->hash = h; s->patternMode = 1;
        }

        /* copy the values of A into the pattern of P (both have sorted row indices within each column) */
        for (j = 0; j < n; j++) {
            mwIndex pk = pjc[j], pend = pjc[j + 1], ak;
            for (ak = pk; ak < pend; ak++) s->pval[ak] = 0.0;
            for (ak = ajc[j]; ak < ajc[j + 1]; ak++) {
                while (pk < pend && pir[pk] < air[ak]) pk++;
                if (pk == pend || pir[pk] != air[ak]) {          /* entry of A outside the pattern of P */
                    if (nlhs > 0) plhs[0] = mxCreateDoubleScalar(-1.0);
                    s->factorised = 0;
                    return;
                }
                s->pval[pk] = aval[ak];
                pk++;
            }
        }

        if (force && s->analysed) { callPardisoRaw(s, -1, NULL, NULL, NULL, NULL, NULL, 1); memset(s->pt, 0, sizeof(s->pt)); s->analysed = 0; }
        if (!s->analysed) {
            callPardisoRaw(s, 11, s->pval, s->pia, s->pja, NULL, NULL, 1);
            s->analysed = 1; reanalysed = 1;
        }
        callPardisoRaw(s, 22, s->pval, s->pia, s->pja, NULL, NULL, 1);
        s->factorised = 1;
        if (nlhs > 0) plhs[0] = mxCreateDoubleScalar((double) reanalysed);
        return;
    }

    if (strcmp(cmd, "solvep") == 0) {
        SolverState *s; i64 nr; double *bcopy, *x;
        requireLoaded();
        if (nrhs < 3) mexErrMsgIdAndTxt("pardisoMex:usage", "x = pardisoMex('solvep', h, b)");
        s = getSolver(prhs[1]);
        if (!s->patternMode || !s->factorised) mexErrMsgIdAndTxt("pardisoMex:state", "Call 'factorp' first.");
        if (mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (i64) mxGetM(prhs[2]) != s->n)
            mexErrMsgIdAndTxt("pardisoMex:rhs", "b must be a real, full, double matrix with as many rows as A.");
        nr = (i64) mxGetN(prhs[2]);
        bcopy = (double *) mxMalloc((size_t) (s->n * nr) * sizeof(double));
        memcpy(bcopy, mxGetDoubles(prhs[2]), (size_t) (s->n * nr) * sizeof(double));
        plhs[0] = mxCreateDoubleMatrix((mwSize) s->n, (mwSize) nr, mxREAL);
        x = mxGetDoubles(plhs[0]);
        callPardisoRaw(s, 33, s->pval, s->pia, s->pja, bcopy, x, nr);
        mxFree(bcopy);
        return;
    }

    mexErrMsgIdAndTxt("pardisoMex:usage", "Unknown command '%s'.", cmd);
}
