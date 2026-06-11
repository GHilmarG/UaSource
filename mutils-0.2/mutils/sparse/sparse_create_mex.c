/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include <libutils/config.h>
#include "opts.h"

#include <libutils/mtypes.h>
#include <libutils/utils.h>
#include <libutils/sorted_list.h>
#include <libutils/memutils.h>
#include <libutils/tictoc.h>

#include <matlab/mexparams.h>
#include <matlab/message_id.h>

/* dof_map uses 1-based numbering */
void sparse_matrix_create(Int matrix_dim, Int n_node_dof,   Int n_elem_dof, 
			  Int n_elems,    Int n_elem_nodes, Int *elems,
			  indexType **Ap_out, dimType **Ai_out, indexType *nnz_out,
			  int symm, int n_row_entries, Int *dof_map)
{
  Int *lists_static = 0;
  Int *n_list_elems_static = 0;

  Int **lists_dynamic = 0;
  Int *n_list_elems_dynamic = 0;
  Int *list_size_dynamic = 0;

  Int i, j, iel;
  Int dofi, dofj;
  Int pos;
  Int *element_dofs = alloca(n_elem_dof*sizeof(Int));

  indexType *Ap = 0;

  /* allocate static data structure first */
  /* if we run out of space, allocate dynamic lists */
  mmalloc(lists_static, sizeof(Int)*n_row_entries*matrix_dim);
  mcalloc(n_list_elems_static, sizeof(Int)*matrix_dim);

  /* allocate storage for dynamic data structure */
  mmalloc(lists_dynamic, sizeof(Int*)*matrix_dim);
  mcalloc(n_list_elems_dynamic, sizeof(Int)*matrix_dim);
  mcalloc(list_size_dynamic, sizeof(Int)*matrix_dim);

  DMESSAGE("static memory usage (MB): %lli\n", DEBUG_MEMORY, (long long)get_total_memory_usage()/1024/1024);

  /* element loop */
  for(iel=0; iel<n_elems; iel++){

    /* enumerate elems dofs */
    for(i=0; i<n_elem_nodes; i++){
      for(j=0; j<n_node_dof; j++){
	dofi = n_node_dof*(elems[iel*n_elem_nodes+i] - ONE_BASED_INDEX)+j;
	if(dof_map){
	  dofi = dof_map[dofi] - ONE_BASED_INDEX;
	}	  
	element_dofs[n_node_dof*i+j] = dofi;
      }
    }

    for(i=0; i<n_elem_dof; i++){

      /* Note the dof traversal order for symmetric matrices */
      for(j=symm*i; j<n_elem_dof; j++){

	dofi = element_dofs[i];
	dofj = element_dofs[j];
	if(symm && dofi>dofj){
	  pos  = dofj;
	  dofj = dofi;
	  dofi = pos;
	}

	/* choose between the static and dynamic data structure */
	if(n_list_elems_static){
	  if(n_list_elems_static[dofi] < n_row_entries){

	    /* Add to static list since there is space  */
	    /* Duplicate entry is not added to the list */
	    sorted_list_add_static_Int(lists_static + dofi*n_row_entries, n_list_elems_static + dofi, dofj);
	    continue;
	  } else {

	    /* locate to see if we have it in the static list already */
	    pos = sorted_list_locate_Int(lists_static + dofi*n_row_entries, n_row_entries, dofj);
	    if(pos<n_row_entries && lists_static[dofi*n_row_entries + pos]==dofj)
	      continue;
	  }
	}

	/* dynamic data structure */
	/* check if a given row has a dynamic list */
	if(!list_size_dynamic[dofi]){
	  sorted_list_create(lists_dynamic+dofi, list_size_dynamic+dofi);
	}
	sorted_list_add(lists_dynamic + dofi, n_list_elems_dynamic + dofi, list_size_dynamic + dofi, dofj);
      }
    }
  }

  DMESSAGE("static and dynamic memory usage (MB): %lli\n", DEBUG_MEMORY, 
	   (long long)get_total_memory_usage()/1024/1024);

  /* allocate output row pointer */
  mmalloc_global(Ap, sizeof(indexType)*(matrix_dim+1));

  /* count number of non-zeros in every matrix row */
  Ap[0] = 0;
  for(i=0; i<matrix_dim; i++) {
    Ap[i+1] = Ap[i] + n_list_elems_static[i] + n_list_elems_dynamic[i];
  }

  /* allocate row array */
  mmalloc_global(*Ai_out, sizeof(dimType)*Ap[matrix_dim]);
  DMESSAGE("maximum memory usage (MB): %lli\n", DEBUG_MEMORY, 
	   (long long)get_total_memory_usage()/1024/1024);

  /* create row entries */

#ifdef SPARSE_CREATE_ALTERNATIVES
  {
    /* This implementation is a standard bucket merging */
    /* Note: this version is slightly slower than the next one. */

    Int optr = 0; /* output data pointer */
    for(i=0; i<matrix_dim; i++){

      Int sptr = 0; /* static data pointer */
      Int dptr = 0; /* dynamic data pointer */

      /* merge sorted lists: static and dynamic */
      for(j=0; j<(Ap[i+1]-Ap[i]); j++){
	if(dptr < n_list_elems_dynamic[i] &&
	   (sptr == n_row_entries || (lists_static[i*n_row_entries + sptr] > lists_dynamic[i][dptr]))){
	  (*Ai_out)[optr] = lists_dynamic[i][dptr];
	  optr++;
	  dptr++;
	} else {
	  (*Ai_out)[optr] = lists_static[i*n_row_entries + sptr];
	  optr++;
	  sptr++;
	}
      }
      if(list_size_dynamic[i]) mfree(lists_dynamic[i], sizeof(Int)*(list_size_dynamic[i]));
    }
  }
#else 
  {
    /* This implementation first copies static list entries */
    /* to the correct row, and then treats the row as a sorted list */
    /* to add the dynamic entries, if any. */

    int ptr;
    dimType *rowptr;
    Int *dynptr;
    dimType n_entries;
    dimType *Ai = *Ai_out;
    for(i=0; i<matrix_dim; i++){

      rowptr = Ai+Ap[i];

      /* copy static entries first */
      n_entries = n_list_elems_static[i];
      for(ptr=0; ptr<n_entries; ptr++){
	rowptr[ptr] = lists_static[i*n_row_entries + ptr];
      }
      
      /* add dynamic entries as to a sorted list */
      if(n_list_elems_dynamic[i]){
	dynptr = lists_dynamic[i];
	for(ptr=0; ptr<n_list_elems_dynamic[i]; ptr++){
	  sorted_list_add_static_dimType(rowptr, &n_entries, dynptr[ptr]);
	}
	mfree(lists_dynamic[i], sizeof(Int)*(list_size_dynamic[i]));
      }
    }
  }
#endif /* SPARSE_CREATE_ALTERNATIVES */
  
  mfree(lists_dynamic, sizeof(Int*)*matrix_dim);
  mfree(n_list_elems_dynamic, sizeof(Int)*(matrix_dim));
  mfree(list_size_dynamic, sizeof(Int)*matrix_dim);

  mfree(lists_static, sizeof(Int)*n_row_entries*matrix_dim);
  mfree(n_list_elems_static, sizeof(Int)*matrix_dim);

  *Ap_out  = Ap;
  *nnz_out = Ap[matrix_dim];
}



void sparse_matrix_create_accum(Int matrix_dim, Int n_node_dof, Int n_elem_dof, 
				Int n_elems, Int n_elem_nodes, Int *elems, Double *Aelems, Int elem_matrices,
				indexType **Ap_out, dimType **Ai_out, Double **Ax_out, indexType *nnz_out, int symm, int n_row_entries,
				Int *dof_map)
{
  Int *lists_static = 0;
  Double *dlists_static = 0;
  Int *n_list_elems_static = 0;

  Int **lists_dynamic = 0;
  Double **dlists_dynamic = 0;
  Int *n_list_elems_dynamic = 0;
  Int *list_size_dynamic = 0;

  indexType *Ap = 0;

  Int i, j, iel;
  Int dofi, dofj;
  Int *element_dofs = alloca(n_elem_dof*sizeof(Int));
  Int pos = 0;
  Double Avalue, *Aelems_ptr = Aelems;

  /* allocate static data structure first */
  /* if we run out of space, allocate dynamic lists */
  mmalloc(lists_static, sizeof(Int)*n_row_entries*matrix_dim);
  mmalloc(dlists_static, sizeof(Double)*n_row_entries*matrix_dim);
  mcalloc(n_list_elems_static, sizeof(Int)*matrix_dim);

  /* allocate storage for dynamic data structure */
  mmalloc(lists_dynamic, sizeof(Int*)*matrix_dim);
  mmalloc(dlists_dynamic, sizeof(Double*)*matrix_dim);
  mcalloc(list_size_dynamic, sizeof(Int)*matrix_dim);
  mcalloc(n_list_elems_dynamic, sizeof(Int)*matrix_dim);

  DMESSAGE("static memory usage (MB): %lli\n", DEBUG_MEMORY, 
	   (long long)get_total_memory_usage()/1024/1024);
  DMESSAGE("distinct elem matrices %d\n", DEBUG_DETAILED, elem_matrices);

  /* element loop */
  for(iel=0; iel<n_elems; iel++){

    /* using the same element matrix for all elements */
    if(!elem_matrices) Aelems_ptr = Aelems;

    /* enumerate elems dofs */
    for(i=0; i<n_elem_nodes; i++){
      for(j=0; j<n_node_dof; j++){
	dofi = n_node_dof*(elems[iel*n_elem_nodes+i] - ONE_BASED_INDEX)+j;
	if(dof_map){
	  dofi = dof_map[dofi] - ONE_BASED_INDEX;
	}
	element_dofs[n_node_dof*i+j] = dofi;
      }
    }

    for(i=0; i<n_elem_dof; i++){

      /* Note the dof traversal order for symmetric matrices */
      for(j=symm*i; j<n_elem_dof; j++){
	
	dofi = element_dofs[i];
	dofj = element_dofs[j];
	if(symm && dofi>dofj){
	  pos  = dofj;
	  dofj = dofi;
	  dofi = pos;
	}

	Avalue = *Aelems_ptr;
	Aelems_ptr++;
	  
	/* choose between the static and dynamic data structure */
	if(n_list_elems_static){
	  
	  /* add to static list since there is space */
	  if(n_list_elems_static[dofi] < n_row_entries){
	    sorted_list_add_static_accum_Int_Double(lists_static + dofi*n_row_entries, n_list_elems_static + dofi, 
						    dofj, dlists_static + dofi*n_row_entries, Avalue);
	    continue;
	  } else {

	    /* locate to see if we have it in the static list already */
	    pos = sorted_list_locate_Int(lists_static + dofi*n_row_entries, n_row_entries, dofj);
	    if(pos<n_row_entries && lists_static[dofi*n_row_entries + pos]==dofj) {
	      dlists_static[dofi*n_row_entries + pos] += Avalue;
	      continue;
	    }
	  }
	} 

	/* dynamic data structure */
	/* check if a given row has a dynamic list */
	if(!list_size_dynamic[dofi]){
	  sorted_list_create_pair(lists_dynamic + dofi, dlists_dynamic + dofi, list_size_dynamic + dofi);
	}
	sorted_list_add_accum(lists_dynamic + dofi, n_list_elems_dynamic + dofi, list_size_dynamic + dofi, dofj,
			      dlists_dynamic + dofi, Avalue);
      }
    }
  }
  DMESSAGE("static and dynamic memory usage (MB): %lli\n", DEBUG_MEMORY,
	   (long long)get_total_memory_usage()/1024/1024);

  /* allocate output row pointer */
  mmalloc_global(Ap, sizeof(indexType)*(matrix_dim+1));
  Ap[0] = 0;

  /* count number of non-zeros */
  for(i=0; i<matrix_dim; i++) {
    Ap[i+1] = Ap[i] + n_list_elems_static[i] + n_list_elems_dynamic[i];
  }

  /* allocate row array */
  mmalloc_global(*Ai_out, sizeof(dimType)*Ap[matrix_dim]);
  mmalloc_global(*Ax_out, sizeof(Double)*Ap[matrix_dim]);
  DMESSAGE("maximum memory usage (MB): %lli\n", DEBUG_MEMORY,
	   (long long)get_total_memory_usage()/1024/1024);

  /* create rows */
#ifdef SPARSE_CREATE_ALTERNATIVES
  {
    Int optr = 0;
    for(i=0; i<matrix_dim; i++){

      Int sptr = 0;
      Int dptr = 0;

      /* merge buckets */
      for(j=0; j<(Ap[i+1]-Ap[i]); j++){
	if(dptr < n_list_elems_dynamic[i] && (sptr == n_row_entries || (lists_static[i*n_row_entries + sptr] > lists_dynamic[i][dptr]))){
	  (*Ai_out)[optr] = lists_dynamic[i][dptr];
	  (*Ax_out)[optr] = dlists_dynamic[i][dptr];
	  optr++;
	  dptr++;
	} else {
	  (*Ai_out)[optr] = lists_static[i*n_row_entries + sptr];
	  (*Ax_out)[optr] = dlists_static[i*n_row_entries + sptr];
	  optr++;
	  sptr++;
	}
      }
      if(list_size_dynamic[i]){
	mfree(lists_dynamic[i], sizeof(Int)*(list_size_dynamic[i]));
	mfree(dlists_dynamic[i], sizeof(Double)*(list_size_dynamic[i]));
      }
    }
  }
#else
  {

    /* This implementation first copies static list entries */
    /* to the correct row, and then treats the row as a sorted list */
    /* to add the dynamic entries, if any. */

    Int ptr;
    dimType *rowptr;
    Double  *rowptr_Ax;
    Int *dynptr;
    Double *ddynptr;
    dimType n_entries;
    dimType *Ai = *Ai_out;
    Double *Ax = *Ax_out;
    for(i=0; i<matrix_dim; i++){

      rowptr = Ai+Ap[i];
      rowptr_Ax = Ax+Ap[i];

      /* copy static entries first */
      n_entries = n_list_elems_static[i];
      for(ptr=0; ptr<n_entries; ptr++){
	rowptr[ptr] = lists_static[i*n_row_entries + ptr];
      }
      memcpy(rowptr_Ax, dlists_static+i*n_row_entries, sizeof(Double)*n_entries);
      
      /* add dynamic entries as to a sorted list */
      if(n_list_elems_dynamic[i]){
	dynptr = lists_dynamic[i];
	ddynptr = dlists_dynamic[i];
	for(ptr=0; ptr<n_list_elems_dynamic[i]; ptr++){
	  sorted_list_add_static_accum_dimType_Double(rowptr, &n_entries, dynptr[ptr], rowptr_Ax, ddynptr[ptr]);
	}
	mfree(lists_dynamic[i], sizeof(Int)*(list_size_dynamic[i]));
	mfree(dlists_dynamic[i], sizeof(Double)*(list_size_dynamic[i]));
      }
    }
  }
#endif
  
  mfree(lists_dynamic, sizeof(Int*)*matrix_dim);
  mfree(dlists_dynamic, sizeof(Double*)*matrix_dim);
  mfree(list_size_dynamic, sizeof(Int)*matrix_dim);
  mfree(n_list_elems_dynamic, sizeof(Int)*(matrix_dim));

  mfree(lists_static, sizeof(Int)*n_row_entries*matrix_dim);
  mfree(dlists_static, sizeof(Double)*n_row_entries*matrix_dim);
  mfree(n_list_elems_static, sizeof(Int)*matrix_dim);

  *Ap_out  = Ap;
  *nnz_out = Ap[matrix_dim];
}


void sparse_map_create(Int n_node_dof, Int n_elem_dof, 
		       Int n_elems, Int n_elem_nodes, Int *elems, Int *dof_map, Int symm,
		       indexType *Ap, dimType *Ai, Int **Map_out, Int *map_size){
  Int  i, j, iel;
  Int *element_dofs = alloca(n_elem_dof*sizeof(Int));
  Int sparse_iter = 0;
  Int sparse_arr_length;
  Int *map;
  Int dofi, dofj, pos;
  indexType *rowptr = 0;
  indexType rowstart = 0;
  indexType nrowent = 0;

  if(symm){
    sparse_arr_length = (Int)(n_elems*0.5*(n_elem_dof)*(n_elem_dof+1));
  } else {
    sparse_arr_length = (Int)(n_elems*(n_elem_dof)*(n_elem_dof));
  }

  mmalloc_global(map, sizeof(Int)*sparse_arr_length);

  for(iel=0; iel<n_elems; iel++){

    for(i=0; i<n_elem_nodes; i++){
      for(j=0; j<n_node_dof; j++){
	dofi = n_node_dof*(elems[iel*n_elem_nodes+i] - ONE_BASED_INDEX)+j;
	if(dof_map){
	  dofi = dof_map[dofi] - ONE_BASED_INDEX;
	}
	element_dofs[n_node_dof*i+j] = dofi;
      }
    }

    for(i=0; i<n_elem_dof; i++){
      for(j=symm*i; j<n_elem_dof; j++){
	  
	dofi = element_dofs[i];
	dofj = element_dofs[j];
	if(symm && dofi>dofj){
	  pos  = dofj;
	  dofj = dofi;
	  dofi = pos;
	}

	rowptr   = Ai + Ap[dofi];
	rowstart = Ap[dofi];
	nrowent  = Ap[dofi+1] - Ap[dofi];

	map[sparse_iter++] = ONE_BASED_INDEX + 
	  (Int)(rowstart + sorted_list_locate_indexType(rowptr, nrowent, dofj));
      }
    }
  }

  *Map_out = map;
  *map_size = sparse_arr_length;
}


void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[])
{
  ti32 n_elem_nodes;
  ti32 n_elems;
  Double *Aelems = NULL;
  ti32 symbolic = 1;
  Int *dof_map = NULL;
  Int *elems = NULL;
  t_opts opts;
  Int distinct_elem_matrices = 1; /* unique element matrix for all elements */
  int arg = 0;
  
  indexType *Ap = NULL;
  dimType *Ai = NULL;
  Double *Ax = NULL;
  indexType matrix_nz;
  Int matrix_dim = 0;
  int argout = 0;

  if (nargin < 1 || nargin > 4) 
    USERERROR("Usage: A = sparse_create(ELEMS, [Aelems=true], [opts], [dof_map=empty])", MUTILS_INVALID_PARAMETER);

  set_debug_mode(0);

  /* ELEMS */
  n_elems = 0;
  n_elem_nodes = 0;
  elems = mex_get_matrix_Int(pargin[arg++], &n_elem_nodes, &n_elems, 
			     "ELEMS", "number of element nodes", 
			     "number of elements", 0);

  /* options */
  arg = 2;
  if(nargin>=arg+1){
    opts = mex2opts(pargin[arg]);
  } else {
    opts = mex2opts(NULL);
  }

  if(opts.n_row_entries == -1){
    switch(n_elem_nodes){
    case  3: opts.n_row_entries = 10; break;
    case  6: opts.n_row_entries = 20; break;
    case  7: opts.n_row_entries = 26; break;
    case  4: opts.n_row_entries = 14; break; /* assume 3d tet, not 2d quad */
    case  9: opts.n_row_entries = 14; break;
    case 10: opts.n_row_entries = 48; break;
    case 15: opts.n_row_entries = 32; break;
    case  8: opts.n_row_entries = 14; break;
    case 27: opts.n_row_entries = 48; break;
    default: opts.n_row_entries = 16; break;
    }
    opts.n_row_entries *= opts.n_node_dof;
  }

  /* element matrix for node dofs */
  arg = 1;
  if(nargin>=arg+1){
    Int m, n;
    char buff[256] = {0};
    SNPRINTF(buff, 255, "Maximum dimensions of Aelem array are nentries X %lli", (long long)MaxInt);
    Aelems = mxGetData(pargin[arg]);
    managed_type_cast(Int, m, mxGetM(pargin[arg]), buff);
    managed_type_cast(Int, n, mxGetN(pargin[arg]), buff);

    symbolic = 0;

    /* decipher Aelems size */
    if(m==1 && n==1){

      /* OK. symbolic matrix */
      symbolic = 1;

    } else if(n==n_elems || m==n_elems){
      if(opts.symmetric){
	if(m*n != ((long)n_elems)*(n_elem_nodes*opts.n_node_dof)*(n_elem_nodes*opts.n_node_dof+1)/2){
	  USERERROR("For symmetric matrices Aelems must be the size of (nnod*ndof)*(nnod*ndof+1)/2 X (1 or nel)", MUTILS_INVALID_PARAMETER);
	}
      } else {
	if(m*n != ((long)n_elems)*(n_elem_nodes*opts.n_node_dof)*(n_elem_nodes*opts.n_node_dof)){
	  USERERROR("For general matrices Aelems must be the size of (nnod*ndof)*(nnod*ndof) X (1 or nel)", MUTILS_INVALID_PARAMETER);
	}	
      }
      /* OK. element matrices for every element separately */
      distinct_elem_matrices = 1;

    } else if(m==1 || n==1){
      if(opts.symmetric){
	/* symmetric sparse matrix and common element matrix */
	if(m*n != (n_elem_nodes*opts.n_node_dof)*(n_elem_nodes*opts.n_node_dof+1)/2){
	  USERERROR("For symmetric matrices Aelems must be the size of (nnod*ndof)*(nnod*ndof+1)/2 X (1 or nel)", MUTILS_INVALID_PARAMETER);
	}
      } else {
	/* general sparse matrix and common element matrix */
	if(m*n != (n_elem_nodes*opts.n_node_dof)*(n_elem_nodes*opts.n_node_dof)){
	  USERERROR("For general matrices Aelems must be the size of (nnod*ndof)*(nnod*ndof) X (1 or nel)", MUTILS_INVALID_PARAMETER);
	}
      }
      /* OK. The same element matrix for all elements */
      distinct_elem_matrices = 0;     
    } else {
      USERERROR("Can not understand size of Aelem.", MUTILS_INVALID_PARAMETER);
    }
  }


  /* create sparse matrix */
  /* find out the matrix dimensions */
  /* TODO add matrix_dim parameter */
  {
    Int i;
    for(i=0; i<n_elems*n_elem_nodes; i++) 
      matrix_dim = matrix_dim<elems[i]?elems[i]:matrix_dim;
    matrix_dim = matrix_dim*opts.n_node_dof;
  }


  /* check if we have a permutation/map */
  arg = 3;
  if(nargin>=arg+1){
    Int m = 1, i;
    dof_map = mex_get_matrix_Int(pargin[arg++], &m, &matrix_dim, 
			     "dof_map", "1", "matrix dimension", 1);
    
    if(dof_map){

      /* dof_map can map nodes/dofs onto each other and reduce matrix_dim */
      m = 0;
      for(i=0; i<matrix_dim; i++) m = MAX(m, dof_map[i]);
      matrix_dim = m;
    }
  }

  /* assemble the matrix, or just create the sparsity structure */
  if(!symbolic){
    sparse_matrix_create_accum(matrix_dim, opts.n_node_dof, opts.n_node_dof*n_elem_nodes, n_elems, n_elem_nodes, 
			       elems, Aelems, distinct_elem_matrices, &Ap, &Ai, &Ax, 
			       &matrix_nz, opts.symmetric, opts.n_row_entries, dof_map);
  } else {
    sparse_matrix_create(matrix_dim, opts.n_node_dof, opts.n_node_dof*n_elem_nodes, n_elems, n_elem_nodes, 
			 elems, &Ap, &Ai, &matrix_nz, opts.symmetric, opts.n_row_entries, dof_map);
  }

  /* return sparse matrix */
  {
    mxArray *A = 0;
    Long i;

    if(!Ax){
      A = mxCreateSparseLogicalMatrix (0, 0, 0);
      mmalloc_global(Ax, matrix_nz*sizeof(mxLogical));
      for(i=0; i<matrix_nz; i++) ((mxLogical*)Ax)[i] = 1;
    } else {
      A = mxCreateSparse (0, 0, 0, mxREAL);
    }

    mpersistent(Ap);
    mpersistent(Ai);
    mpersistent(Ax);

    mxSetM(A, matrix_dim);
    mxSetN(A, matrix_dim);
    mxSetNzmax(A, matrix_nz);
    mxSetJc(A, Ap);
    mxSetIr(A, Ai);
    mxSetPr(A, Ax);

    pargout[argout] = A;

    if(symbolic)
      dec_memory_usage(matrix_nz*(sizeof(mxLogical)+sizeof(dimType)));
    else
      dec_memory_usage(matrix_nz*(sizeof(Double)+sizeof(dimType)));
    dec_memory_usage((matrix_dim+1)*sizeof(indexType));
  }

  /* create map */
  if(opts.gen_map && nargout>1){
    Int *map = NULL, map_size = 0;
    Int n_elem_dof = n_elem_nodes*opts.n_node_dof;
    Int size;
    sparse_map_create(opts.n_node_dof, opts.n_node_dof*n_elem_nodes, n_elems, n_elem_nodes,
  		      elems, dof_map, opts.symmetric, Ap, Ai, &map, &map_size);
    pargout[1] = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
    mxSetN(pargout[1], 1);
    mxSetM(pargout[1], map_size);
    mxSetData(pargout[1], map);
    mpersistent(map);

    /* 'free' the allocated memory that we release to MATLAB */
    /* not our responsibility anymore */
    if(opts.symmetric){
      size = (Int)(n_elems*0.5*(n_elem_dof)*(n_elem_dof+1));
    } else {
      size = (Int)(n_elems*(n_elem_dof)*(n_elem_dof));
    }
    dec_memory_usage(sizeof(Int)*size);
  }
}
