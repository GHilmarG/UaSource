#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <sys/mman.h>
#include <strings.h>
#include <string.h>

/* typedef __int32_t  mwSize; */
/* typedef __uint32_t mwIndex; */

#define MAX(a, b) (a)>(b)?(a):(b)
#define MIN(a, b) (a)<(b)?(a):(b)

typedef __int32_t mint;
typedef __int32_t muint;
typedef __int32_t mint32;

#define mmalloc(v, size) \
  {			 \
    v = malloc(size);	 \
  }			 \

#define mcalloc(v, size) \
  {			 \
    v=malloc(size);	 \
    bzero(v, size);	 \
  }

#define mfree(v, size)	 \
  {			 \
    free(v);		 \
  }

#define mrealloc(v, size)	 \
  {				 \
    v = realloc(v, size);	 \
  }

struct timeval tb, te;
long long flops;
void tic()
{
  gettimeofday(&tb, NULL);
  flops=0;
  fflush(stdout);
}

void toc()
{
  long s,u;
  double tt;
  gettimeofday(&te, NULL);
  s=te.tv_sec-tb.tv_sec;
  u=te.tv_usec-tb.tv_usec;
  tt=((double)s)*1000000+u;
  fprintf(stderr, "time: %li.%.6li\n", (s*1000000+u)/1000000, (s*1000000+u)%1000000);
}


mint triplet2crs(muint n_row, muint n_col, muint nz, muint n_elem_nodes, muint n_elems, 
		 mint32 Ti[], mint32 Tj[], mint32 elems[], double Ax[], muint symm, muint n_node_dof, 
		 muint **Ap_out, muint **Ai_out, double **Ax_out, muint *matrix_nz)
{
  muint i, j, k, iel, temp, l, u;
  muint *buckets_r, *buckets_rz, *buckets_rzz;
  muint *buckets_c, *buckets_cz;
  muint *ncolent;
  muint *temp_c;
  muint *rowform;
  muint *colform;
  double *Ax_rowform;
  double *Ax_colform;
  muint dest, row, col;
  muint dof1, dof2, dofoffset;
  double val;
  muint element_dofs[n_elem_nodes*n_node_dof];

  /* column buckets */
  mcalloc(buckets_cz,  sizeof(muint)*(n_col+1));
  buckets_c   = buckets_cz + 1;
  
  /* row buckets */
  mcalloc(buckets_rzz, sizeof(muint)*(n_row+2));
  buckets_rz   = buckets_rzz + 1;
  buckets_r    = buckets_rzz + 2;

  mmalloc(temp_c,       sizeof(muint)*(n_col+1));
  mcalloc(ncolent,      sizeof(muint)*(n_col+1));
  
  tic();

  /* calculate row buckets - number of entries */
  /* in each row */
  if(elems){

    /* based on the element structure */
    for(iel=0; iel<n_elems; iel++){

      /* find element dofs */
      for(i=0; i<n_elem_nodes; i++){
	for(j=0; j<n_node_dof; j++){
	  element_dofs[n_node_dof*i+j] =  n_node_dof*(elems[iel*n_elem_nodes+i]-1)+j;
	}
      }

      for(j=0; j<n_elem_nodes*n_node_dof; j++){
	dof1 = element_dofs[j];
	for(k=symm*j; k<n_elem_nodes*n_node_dof; k++){
	  dof2 = element_dofs[k];
	  if(dof1>dof2 || !symm){
	    buckets_r[dof1]++;
	    buckets_c[dof2]++;
	  } else {
	    buckets_r[dof2]++;
	    buckets_c[dof1]++;
	  }
	}
      }
    }
  } else {
  
    /* based on the triplet index arrays */
    for(i=0; i<nz; i++) buckets_r[Ti[i]-1]++;
    for(i=0; i<nz; i++) buckets_c[Tj[i]-1]++;
  }

  /* calculate bucket pointers in a contiguous */
  /* array (sum up bucket sizes) */
  for(i=1; i<n_row; i++) buckets_r[i] += buckets_r[i-1];
  for(i=0; i<n_col; i++) buckets_c[i] += buckets_c[i-1];
  toc();

  /* allocate memory for sparse structure, including duplicates */
  nz = buckets_rz[n_row];
  mmalloc(rowform, nz*sizeof(muint));
  mmalloc(colform, nz*sizeof(muint));
  if(Ax) {
    mmalloc(Ax_rowform, nz*sizeof(double));
    mmalloc(Ax_colform, nz*sizeof(double));
  }
  
  /*
   * copy column value of every element to proper buckets
   * ~sort by row
   */
  tic();
  if(elems){
    
    int Aelem_entries = (n_elem_nodes*n_node_dof)*(n_elem_nodes*n_node_dof+symm)/(2*symm);

   /* based on the element structure */
    for(iel=0; iel<n_elems; iel++){

      /* find element dofs */
      for(i=0; i<n_elem_nodes; i++){
	for(j=0; j<n_node_dof; j++){
	  element_dofs[n_node_dof*i+j] =  n_node_dof*(elems[iel*n_elem_nodes+i]-1)+j;
	}
      }

      dofoffset = 0;
      for(j=0; j<n_elem_nodes*n_node_dof; j++){
	dof1 = element_dofs[j];
	for(k=symm*j; k<n_elem_nodes*n_node_dof; k++){
	  dof2 = element_dofs[k];
	  if(dof1>dof2 || !symm){
	    dest = buckets_rz[dof1]++;
	    rowform[dest] = dof2;
	  } else {
	    dest = buckets_rz[dof2]++;
	    rowform[dest] = dof1;
	  }
	  if(Ax) Ax_rowform[dest] = Ax[iel*Aelem_entries + dofoffset + k-j];
	}
	dofoffset += n_elem_nodes*n_node_dof - j;
      }
    }
  } else {

    /* based on the triplet index arrays */
    for(i=0; i<nz; i++) {
      dest = buckets_rz[Ti[i]-1]++;
      rowform[dest] = Tj[i]-1;
      if(Ax) Ax_rowform[dest] = Ax[i];
    }
  }
  toc();

  /* remove duplicates and transpose - sort by columns */
  tic();
  for(i=0; i<n_col; i++) temp_c[i] = (muint)(-1);
  for(i=0; i<n_row; i++){

    row = i;
    l   = buckets_rzz[i];
    u   = buckets_rzz[i+1];

    for(j=l; j<u; j++) {

      col = rowform[j];
      if(Ax) val = Ax_rowform[j];

      if(temp_c[col] != row){

	/* this row entry not yet present in the column */
	temp_c[col] = row;
	dest = buckets_cz[col]++;
	colform[dest] = row;
	if(Ax) Ax_colform[dest] = val;
	ncolent[col]++;
      } else {

	/* row already added, only update Ax */
	if(Ax) Ax_colform[buckets_cz[col]-1] += val;
      }
    }
  }
  toc();

  /* compact the sparse storage */
  tic();
  dest           = ncolent[0];
  ncolent[0]     = 0;
  buckets_cz[0] -= dest;
  for(i=1; i<=n_col; i++){
    temp         = ncolent[i];
    ncolent[i]   = ncolent[i-1]+dest;
    dest         = temp;
    buckets_cz[i] -= dest;    
    memmove(colform + ncolent[i], colform + buckets_cz[i], sizeof(muint)*dest);
    if(Ax) memmove(Ax_colform + ncolent[i], Ax_colform + buckets_cz[i], sizeof(double)*dest);
  }
  toc();
  
  mfree(rowform,   nz*sizeof(muint));
  mfree(buckets_cz,   sizeof(muint)*(n_col+1));
  mfree(temp_c,       sizeof(muint)*(n_col+1));
  mfree(buckets_rzz,  sizeof(muint)*(n_row+2));

  mrealloc(colform, sizeof(muint)*ncolent[n_col]);
  if(Ax){
    mrealloc(Ax_colform, sizeof(double)*ncolent[n_col]);
    mfree(Ax_rowform,   nz*sizeof(double));
    *Ax_out = Ax_colform;
  }

  *Ap_out = ncolent;
  *Ai_out = colform;
  *matrix_nz = ncolent[n_col];

  return 0;
}


void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin [ ])
{
  muint n_elem_nodes = 0;
  muint n_elems = 0;
  muint n_node_dof = 1;
  muint *elems = NULL;
  muint *Ai    = NULL;
  muint *Aj    = NULL;
  double *Ax   = NULL;
  muint arg    = 0;
  muint matrix_nz = 0;
  muint m = 0, n = 0, symbolic=0, symm=0;
  muint Am = 0, An = 0;
  muint n_row = 0;
  muint n_col = 0;
  muint i;

  /* A = mk_sparse5(ELEMS,  [Aelem=true], [n_node_dof=1], [symm=0], [matrix_dim])  */
  /* A = mk_sparse5(Ai, Aj, [Aelem=true], [n_row], [n_col]) */

  if(nargin<1){
    mexErrMsgTxt ("Usage: A = mk_sparse5(ELEMS, [Aelem], [dim], [symm])");
  }

  /* analyze the parameters */
  {
    m = mxGetM(pargin[arg]);
    n = mxGetN(pargin[arg]);
    mxClassID cid = mxGetClassID(pargin[arg]);
    if(!(mxIsUint32(pargin[arg]) || mxIsInt32(pargin[arg]))) 
      mexErrMsgTxt("ELEMS and Ai/Aj must be of type (u)int32");

    /* if a vector, this should be Ai. */
    /* we must also have Aj.           */
    if(m==1 || n==1){
      
      /* Ai Aj syntax */
      if(nargin<2){
	mexErrMsgTxt("Ai/Aj index arrays must be both given.");
      }
      matrix_nz = MAX(m, n);
      Ai    = mxGetData(pargin[arg]);
      arg++;

      /* check Aj */
      m = MAX(mxGetM(pargin[arg]), mxGetN(pargin[arg]));
      n = MIN(mxGetM(pargin[arg]), mxGetN(pargin[arg]));
      mxClassID cid = mxGetClassID(pargin[arg]);
      if(!(mxIsUint32(pargin[arg]) || mxIsInt32(pargin[arg]))) 
	mexErrMsgTxt("ELEMS and Ai/Aj must be of type (u)int32");

      if(n!=1 || m!=matrix_nz)
	mexErrMsgTxt("Ai/Aj dimensions do not match");

      Aj    = mxGetData(pargin[arg]);
      arg++;

      /* next we have Aelem */
      if(nargin>arg){
	Am = mxGetM(pargin[arg]);
	An = mxGetN(pargin[arg]);
	mxClassID cid = mxGetClassID(pargin[arg]);
	if(!(mxIsDouble(pargin[arg]) || mxIsLogical(pargin[arg])))
	  mexErrMsgTxt("Aelem must be of type double or 'true'");
	if(Am*An) Ax = mxGetData(pargin[arg]);
	arg++;
      }

      /* n_row */
      if(nargin>arg){
	m = mxGetM(pargin[arg])*mxGetN(pargin[arg]);
	if(m!=1)
	  mexErrMsgTxt("n_row must be a scalar value");
	if(!(mxIsDouble(pargin[arg])))
	  mexErrMsgTxt("n_row must be of type double");
	n_row = mxGetPr(pargin[arg])[0];
	arg++;
      }

      /* n_col */
      if(nargin>arg){
	m = mxGetM(pargin[arg])*mxGetN(pargin[arg]);
	if(m!=1)
	  mexErrMsgTxt("n_col must be a scalar value");
	if(!(mxIsDouble(pargin[arg])))
	  mexErrMsgTxt("n_col must be of type double");
	n_col = mxGetPr(pargin[arg])[0];
	arg++;
      }

      /* decipher Aelem dimensions if we have it */
      n = An; m = Am;
      if(n*m){
	if(n==matrix_nz || m==matrix_nz){
	  /* OK. Aelem has size of Ai/j */
	} else if(m==1 || n==1){
	  /* OK. Aelem is a scalar.  */	  
	  symbolic = 1;
	  Ax = NULL;
	} else {
	  mexErrMsgTxt("Aelem must be the same size as Ai and Aj arrays, \
or a scalar to create a symbolic sparse matrix.");
	}
      } else symbolic = 1;

      if(!(n_row && n_col)){
	for(i=0; i<matrix_nz; i++){
	  n_row = MAX(n_row, Ai[i]);
	  n_col = MAX(n_row, Aj[i]);
	}
      }
    } else {

      /* ELEMS syntax */
      n_elem_nodes = m;
      n_elems      = n;
      elems        = mxGetData(pargin[arg]);
      arg++;

      /* next we have Aelem */
      if(nargin>arg){
	Am = mxGetM(pargin[arg]);
	An = mxGetN(pargin[arg]);
	mxClassID cid = mxGetClassID(pargin[arg]);
	if(!(mxIsDouble(pargin[arg]) || mxIsLogical(pargin[arg])))
	  mexErrMsgTxt("Aelem must be of type double or 'true'");
	if(Am*An) Ax = mxGetData(pargin[arg]);
	arg++;
      }

      /* n_node_dof */
      if(nargin>arg){
	m = mxGetM(pargin[arg])*mxGetN(pargin[arg]);
	if(m!=1)
	  mexErrMsgTxt("n_node_dof must be a scalar value");
	if(!(mxIsDouble(pargin[arg])))
	  mexErrMsgTxt("n_node_dof must be of type double");
	n_node_dof = mxGetPr(pargin[arg])[0];
	arg++;
      }

      /* symm */
      if(nargin>arg){
	m = mxGetM(pargin[arg])*mxGetN(pargin[arg]);
	if(m!=1)
	  mexErrMsgTxt("symm must be a scalar value of 0 or 1");
	if(!(mxIsDouble(pargin[arg])))
	  mexErrMsgTxt("symm must be of type double");
	symm = mxGetPr(pargin[arg])[0];
	arg++;
      }

      /* matrix_dim */
      if(nargin>arg){
	m = mxGetM(pargin[arg])*mxGetN(pargin[arg]);
	if(m!=1)
	  mexErrMsgTxt("marix_dim must be a scalar value");
	if(!(mxIsDouble(pargin[arg])))
	  mexErrMsgTxt("marix_dim must be of type double");
	n_row = n_col = mxGetPr(pargin[arg])[0];
	arg++;
      }

      /* decipher Aelem dimensions if we have it */
      n = An; m = Am;
      if(n*m){
	if(n==n_elems || m==n_elems){
	  if(m==1 && n==1){
	    /* OK. symbolic matrix */
	    symbolic = 1;
	    Ax = 0;
	  } else {
	    if(symm){
	      if(m*n != ((long)n_elems)*(n_elem_nodes*n_node_dof)*(n_elem_nodes*n_node_dof+1)/2){
		mexErrMsgTxt("For symmetric matrices Aelems must be the size of \
(nnod*ndof)*(nnod*ndof+1)/2 X (1 or nel)");
	      }
	    } else {
	      if(m*n != ((long)n_elems)*(n_elem_nodes*n_node_dof)*(n_elem_nodes*n_node_dof)){
		mexErrMsgTxt("For general matrices Aelems must be the size of \
(nnod*ndof)*(nnod*ndof) X (1 or nel)");
	      }	
	    }
	    /* OK. element matrices for every element separately */
	  }
	} else if(m==1 || n==1){
	  if(m==1 && n==1){
	    /* OK. symbolic matrix */
	    symbolic = 1;
	    Ax = 0;
	  } else {
	    if(symm){
	      /* symmetric sparse matrix and common element matrix */
	      if(m*n != (n_elem_nodes*n_node_dof)*(n_elem_nodes*n_node_dof+1)/2){
		mexErrMsgTxt("For symmetric matrices Aelems must be the size of \
(nnod*ndof)*(nnod*ndof+1)/2 X (1 or nel)");
	      }
	    } else {
	      /* general sparse matrix and common element matrix */
	      if(m*n != (n_elem_nodes*n_node_dof)*(n_elem_nodes*n_node_dof)){
		mexErrMsgTxt("For general matrices Aelems must be the size of \
(nnod*ndof)*(nnod*ndof) X (1 or nel)");
	      }
	    }
	  }
	} else {
	  mexErrMsgTxt("Can not understand size of Aelem.");
	}
      } else symbolic = 1;

      if(!(n_row && n_col)){
	for(i=0; i<n_elems*n_elem_nodes; i++){
	  n_row = MAX(n_row, elems[i]);
	}
	n_row *= n_node_dof;
	n_col = n_row;
      }
    }
  }

  /* read parameters */
  fprintf(stderr, "given n_row %d, n_col %d\n", n_row, n_col);
  fprintf(stderr, "symbolic %d, n_node_dof %d, symm %d, matrix_nz %d, Ax %x, Ai %x, Aj %x, elems %x\n", 
	  symbolic, n_node_dof, symm, matrix_nz, Ax, Ai, Aj, elems);

  muint *Ap_out, *Ai_out;
  double *Ax_out;

  //return;
  triplet2crs(n_row, n_col, matrix_nz, n_elem_nodes, n_elems, Ai, Aj, elems, Ax, 
	      symm, n_node_dof, &Ap_out, &Ai_out, &Ax_out, &matrix_nz);

  if(nargout>0) {
    mxArray *A;
    muint i;
    tic();
    if(symbolic){
      mmalloc(Ax_out, matrix_nz*sizeof(double));
      for(i=0; i<matrix_nz; i++) Ax_out[i] = 1.0;
    }
    toc();

    A = mxCreateSparse (0, 0, 0, mxREAL);
    mxSetM(A, n_row);
    mxSetN(A, n_col);
    mxSetNzmax(A, matrix_nz);

    tic();
    if(sizeof(mwIndex)==8){

      mwIndex *Ap=0;
      mwIndex *Ai=0;

      mmalloc(Ap, sizeof(mwIndex)*(n_row+1));
      for(i=0; i<=n_row; i++) Ap[i] = Ap_out[i];
      mfree(Ap_out, sizeof(muint)*(n_row+1));

      mmalloc(Ai, sizeof(mwIndex)*(matrix_nz));
      for(i=0; i<matrix_nz; i++) Ai[i] = Ai_out[i];
      mfree(Ai_out, sizeof(muint)*(matrix_nz));

      mxSetJc(A, Ap);
      mxSetIr(A, Ai);
      mexMakeMemoryPersistent(Ap);
      mexMakeMemoryPersistent(Ai);
    } else {
      mxSetJc(A, (mwIndex*)Ap_out);
      mxSetIr(A, (mwIndex*)Ai_out);
      mexMakeMemoryPersistent(Ap_out);
      mexMakeMemoryPersistent(Ai_out);
    }
    toc();
    mxSetPr(A, Ax_out);
    mexMakeMemoryPersistent(Ax_out);
    pargout[0] = A;
  } else {
    mfree(Ap_out, 0);
    mfree(Ai_out, 0);
  }
}
