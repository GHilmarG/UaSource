/* SparseFromCSC.c
 *
 *   K = SparseFromCSC(jc, ir, pr, m)
 *
 * Creates an m x n sparse matrix directly from its compressed-column (CSC) arrays, without sorting:
 *   jc : uint64, n+1 zero-based column pointers (jc(1)=0, jc(end)=nnz)
 *   ir : uint64, nnz zero-based row indices, sorted within each column
 *   pr : double, nnz values
 *   m  : number of rows
 * The arrays are taken as they are (they come from the cached assembly pattern, see uvhCSCPattern.m); only their sizes and the first
 * and last column pointer are checked. Explicit zeros in pr are kept as stored entries. The arrays are allocated without zero-filling
 * and copied in parallel (OpenMP, if compiled with OpenMP; otherwise serially). Built by BuildUaCompiled.m.
 * (25 Oct 2026)
 */
#include "mex.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void pcopy(void *dst, const void *src, size_t nbytes)
{
    const size_t chunk = 1u << 20;   /* 1 MB blocks */
    const long long nb = (long long) ((nbytes + chunk - 1) / chunk);
    long long b;
    #pragma omp parallel for schedule(static)
    for (b = 0; b < nb; b++) {
        size_t off = (size_t) b * chunk;
        size_t len = (off + chunk <= nbytes) ? chunk : nbytes - off;
        memcpy((char *) dst + off, (const char *) src + off, len);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4) mexErrMsgIdAndTxt("SparseFromCSC:nrhs", "Four inputs required: jc, ir, pr, m.");
    if (!mxIsUint64(prhs[0]) || !mxIsUint64(prhs[1])) mexErrMsgIdAndTxt("SparseFromCSC:class", "jc and ir must be uint64.");
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])) mexErrMsgIdAndTxt("SparseFromCSC:class", "pr must be a real full double array.");
    if (sizeof(mwIndex) != 8) mexErrMsgIdAndTxt("SparseFromCSC:mwIndex", "mwIndex is not 64 bit.");
    mwSize n  = mxGetNumberOfElements(prhs[0]) - 1;
    mwSize nz = mxGetNumberOfElements(prhs[2]);
    if (mxGetNumberOfElements(prhs[1]) != nz) mexErrMsgIdAndTxt("SparseFromCSC:size", "ir and pr must have the same number of elements.");
    mwSize m  = (mwSize) mxGetScalar(prhs[3]);
    const mwIndex *jcIn = (const mwIndex *) mxGetData(prhs[0]);
    if (jcIn[0] != 0 || jcIn[n] != nz) mexErrMsgIdAndTxt("SparseFromCSC:jc", "jc must start with 0 and end with nnz.");
    mxArray *K = mxCreateSparse(m, n, 1, mxREAL);   /* minimal allocation, replaced below */
    mwIndex *jc = (mwIndex *) mxMalloc((n + 1) * sizeof(mwIndex));
    mwIndex *ir = (mwIndex *) mxMalloc((nz > 0 ? nz : 1) * sizeof(mwIndex));
    double  *pr = (double  *) mxMalloc((nz > 0 ? nz : 1) * sizeof(double));
    pcopy(jc, jcIn, (n + 1) * sizeof(mwIndex));
    pcopy(ir, mxGetData(prhs[1]), nz * sizeof(mwIndex));
    pcopy(pr, mxGetDoubles(prhs[2]), nz * sizeof(double));
    mxFree(mxGetJc(K)); mxFree(mxGetIr(K)); mxFree(mxGetDoubles(K));
    mxSetJc(K, jc); mxSetIr(K, ir); mxSetDoubles(K, pr); mxSetNzmax(K, nz > 0 ? nz : 1);
    plhs[0] = K;
}
