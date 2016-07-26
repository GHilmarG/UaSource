/*
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include <libutils/config.h>
#include <libutils/utils.h>
#include <libutils/parallel.h>
#include <matlab/mesh.h>
#include <matlab/mexparams.h>


void interpolate_in_triangles_reference(t_mesh mesh, double *values, double *markers, Int *element_id, Int n_markers, double *values_markers)
{
  Int i;
  Int elid;
  double area;
  const double *a, *b, *c;
  double eta1, eta2, eta3, eta123;
  double N[7];
  Int n = 0, node;

  for(i=0; i<n_markers; i++){

    elid = element_id[i]-1;

    /* validate the element map */
    if(elid<0 || elid>=mesh.n_elems){
      ERROR("invalid element_id: %d", elid);
    }

    a = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+0]-1);
    b = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+1]-1);
    c = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+2]-1);

    area =
      b[0]*c[1] - c[0]*b[1] +
      c[0]*a[1] - a[0]*c[1] +
      a[0]*b[1] - b[0]*a[1] ;

    eta1 = 
      b[0]*c[1]-c[0]*b[1]+
      c[0]*markers[2*i+1]-markers[2*i+0]*c[1]+
      markers[2*i+0]*b[1]-b[0]*markers[2*i+1];

    eta2 = 
      c[0]*a[1]-a[0]*c[1]+
      a[0]*markers[2*i+1]-markers[2*i+0]*a[1]+
      markers[2*i+0]*c[1]-c[0]*markers[2*i+1];

    eta1 /= area;
    eta2 /= area;
    eta3 = 1-eta1-eta2;

    values_markers[2*i+0] = eta1;
    values_markers[2*i+1] = eta2;
    
    eta123 =  eta1*eta2*eta3;
    N[0] = eta1*(2*eta1-1) + 3*eta123;
    N[1] = eta2*(2*eta2-1) + 3*eta123;
    N[2] = eta3*(2*eta3-1) + 3*eta123;
    N[3] = 4*eta2*eta3 - 12*eta123;
    N[4] = 4*eta1*eta3 - 12*eta123;
    N[5] = 4*eta1*eta2 - 12*eta123;
    N[6] =               27*eta123;
    
    for(n=0; n<mesh.n_elem_nodes; n++){
      node = mesh.elems[elid*mesh.n_elem_nodes + n]-1;
      values_markers[2*i+0] += values[node*2+0]*N[n];
      values_markers[2*i+1] += values[node*2+1]*N[n];
    }
  }
}


#ifndef __SSE2__
#pragma message ("Using non-SSE interpolation.")
void interpolate_in_triangles(t_mesh mesh, double *values, double *markers, Int *element_id, Int n_markers, double *values_markers)
{
  Int i;
  Int elid;
  double area;
  const double *a, *b, *c;
  double eta1, eta2, eta3, eta123;
  double N;
  Int node;
  register double ox, oy;

  for(i=0; i<n_markers; i++){

    elid = element_id[i]-1;

    /* validate the element map */
    if(elid<0 || elid>=mesh.n_elems){
      ERROR("invalid element_id: %d", elid);
      break; /* needed - performance */
    }

    a = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+0]-1);
    b = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+1]-1);
    c = mesh.nodes + 2*(mesh.elems[elid*mesh.n_elem_nodes+2]-1);

    area =
      b[0]*c[1] - c[0]*b[1] +
      c[0]*a[1] - a[0]*c[1] +
      a[0]*b[1] - b[0]*a[1] ;

    eta1 =
      b[0]*c[1]-c[0]*b[1]+
      c[0]*markers[2*i+1]-markers[2*i+0]*c[1]+
      markers[2*i+0]*b[1]-b[0]*markers[2*i+1];

    eta2 =
      c[0]*a[1]-a[0]*c[1]+
      a[0]*markers[2*i+1]-markers[2*i+0]*a[1]+
      markers[2*i+0]*c[1]-c[0]*markers[2*i+1];

    area  = 1.0/area;
    eta1 *= area;
    eta2 *= area;
    eta3  = 1-eta1-eta2;
    eta123 =  eta1*eta2*eta3;

    N = eta1*(2*eta1-1) + 3*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 0]-1;
    ox  = values[node*2+0]*N;
    oy  = values[node*2+1]*N;

    N = eta2*(2*eta2-1) + 3*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 1]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    N = eta3*(2*eta3-1) + 3*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 2]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    N = 4*eta2*eta3 - 12*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 3]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    N = 4*eta1*eta3 - 12*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 4]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    N = 4*eta1*eta2 - 12*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 5]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    N =               27*eta123;
    node = mesh.elems[elid*mesh.n_elem_nodes + 6]-1;
    ox += values[node*2+0]*N;
    oy += values[node*2+1]*N;

    values_markers[2*i+0] = ox;
    values_markers[2*i+1] = oy;
  }
}

#else
#include <emmintrin.h>
#ifdef DEBUG
#pragma message ("Using SSE interpolation.")
#endif
void interpolate_in_triangles(t_mesh mesh, double *values, double *markers, Int *element_id, Int n_markers, double *values_markers)
{

  /* use default/environment defined number of threads */
  parallel_set_num_threads(0);

#ifdef USE_OPENMP
#pragma omp parallel
#endif

  {
    Int i, elid;
    Int marker_start, marker_end;
    double area;
    double eta1, eta2, eta3, eta123;
    Int *nodes, *nodesp1, node;
    __m128d v_val, v_shp, v_res;
    __m128d v_a , v_b , v_c , v_m;
    __m128d v_ar, v_br, v_cr, v_mr;
    __m128d v_bccb, v_caac;
    Int thrid;
    Int nthr;
    Int blk_size;

    parallel_get_info(&thrid, &nthr);

    blk_size     = n_markers/nthr+1;
    marker_start = thrid*blk_size;
    marker_end   = (thrid+1)*blk_size;
    marker_end   = MIN(n_markers, marker_end);

    /* prefetch coordinates for next iteration */
    elid  = element_id[marker_start]-1;

    /* validate the element map */
    if(elid<0 || elid>=mesh.n_elems){
      ERROR("invalid element_id: %d", elid);
    }

    nodesp1 = mesh.elems+elid*mesh.n_elem_nodes;
    v_a = _mm_load_pd(mesh.nodes + 2*(nodesp1[0]-1));
    v_b = _mm_load_pd(mesh.nodes + 2*(nodesp1[1]-1));
    v_c = _mm_load_pd(mesh.nodes + 2*(nodesp1[2]-1));
    v_m = _mm_load_pd(markers + 2*marker_start);

    for(i=marker_start; i<marker_end; i++){

      /* reversed order */
      v_ar = _mm_shuffle_pd(v_a, v_a, _MM_SHUFFLE2 (0,1));
      v_br = _mm_shuffle_pd(v_b, v_b, _MM_SHUFFLE2 (0,1));
      v_cr = _mm_shuffle_pd(v_c, v_c, _MM_SHUFFLE2 (0,1));
      v_mr = _mm_shuffle_pd(v_m, v_m, _MM_SHUFFLE2 (0,1));

      /* area */
      v_caac = _mm_mul_pd(v_c, v_ar);
      v_bccb = _mm_mul_pd(v_b, v_cr);
      v_res  = _mm_mul_pd(v_a, v_br);
      v_res  = _mm_add_pd(v_res, v_caac);
      v_res  = _mm_add_pd(v_res, v_bccb);
      area   = _mm_cvtsd_f64 (_mm_sub_pd(v_res, _mm_shuffle_pd(v_res, v_res, _MM_SHUFFLE2 (0,1))));

      /* SSE3 version not faster */
      /* area   = _mm_cvtsd_f64 (_mm_hsub_pd(v_res, v_res)); */

      /* eta1 */
      v_res = _mm_add_pd(v_bccb, _mm_mul_pd(v_c, v_mr));
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_m, v_br));
      eta1  = _mm_cvtsd_f64 (_mm_sub_pd(v_res, _mm_shuffle_pd(v_res, v_res, _MM_SHUFFLE2 (0,1))));

      /* eta2 */
      v_res = _mm_add_pd(v_caac, _mm_mul_pd(v_a, v_mr));
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_m, v_cr));
      eta2  = _mm_cvtsd_f64 (_mm_sub_pd(v_res, _mm_shuffle_pd(v_res, v_res, _MM_SHUFFLE2 (0,1))));

      area   = 1.0/area;
      eta1  *= area;
      eta2  *= area;
      eta3   = 1-eta1-eta2;
      eta123 = eta1*eta2*eta3;

      /* prefetch coordinates for next iteration */
      nodes = nodesp1;
      if(i<n_markers-1){

	/* save cache on marker coordinates */
	_mm_prefetch((char*)markers+2*(i+1), _MM_HINT_NTA);
  
	elid  = element_id[i+1]-1;

	/* validate the element map */
	if(elid<0 || elid>=mesh.n_elems){
	  ERROR("invalid element_id: %d", elid);
	  break; /* needed - performance */
	}

	nodesp1 = mesh.elems+elid*mesh.n_elem_nodes;
	v_a = _mm_load_pd(mesh.nodes + 2*(nodesp1[0]-1));
	v_b = _mm_load_pd(mesh.nodes + 2*(nodesp1[1]-1));
	v_c = _mm_load_pd(mesh.nodes + 2*(nodesp1[2]-1));
	v_m = _mm_load_pd(markers+2*(i+1));
      }

      /* multiply values at nodes by nodal shape functions in markers */
      node  = nodes[0]-1;
      v_shp = _mm_set1_pd(eta1*(2*eta1-1) + 3*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_mul_pd(v_val, v_shp);

      node  = nodes[1]-1;
      v_shp = _mm_set1_pd(eta2*(2*eta2-1) + 3*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      node  = nodes[2]-1;
      v_shp = _mm_set1_pd(eta3*(2*eta3-1) + 3*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      node  = nodes[3]-1;
      v_shp = _mm_set1_pd(4*eta2*eta3 - 12*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      node  = nodes[4]-1;
      v_shp = _mm_set1_pd(4*eta1*eta3 - 12*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      node  = nodes[5]-1;
      v_shp = _mm_set1_pd(4*eta1*eta2 - 12*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      node  = nodes[6]-1;
      v_shp = _mm_set1_pd(27*eta123);
      v_val = _mm_load_pd(values+node*2);
      v_res = _mm_add_pd(v_res, _mm_mul_pd(v_val, v_shp));

      /* stream the result */
      _mm_stream_pd(values_markers+2*i, v_res);
    }
  }
}
#endif


void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[])
{
  Int arg = 0;
  Int n_dim, temp;

  Int n_values;
  Double *values;
  Double *values_markers;

  Int n_markers;
  Double *markers;

  Int *element_id;
  mxArray *outp;
  t_mesh mesh;

  if (nargin < 4){
    mexErrMsgTxt("Usage: V_MARKERS = einterp(MESH, V, MARKERS, element_id)");
  }

  /* analyze arguments */
  mesh = mex2mesh(pargin[arg++]);

  /* TODO: for now */
  if(mesh.n_elem_nodes!=7)
    mexErrMsgTxt("Only works for 7-node triangles for now");

  n_values = mesh.n_nodes;
  n_dim    = 2;
  values   = mex_get_matrix_Double(pargin[arg++], &n_dim, &n_values, "V", "2", "number of mesh nodes", 0);

  n_markers  = 0;
  markers    = mex_get_matrix_Double(pargin[arg++], &n_dim, &n_markers, "MARKERS", "2", "number of markers", 0);

  temp = 1;
  element_id = mex_get_matrix_Int(pargin[arg++], &temp, &n_markers, "ELEMENT_MAP", "1", "number of markers", 0);

  if(nargout==0) return;

  /* create output array */
  outp = mex_set_matrix_Double(NULL, 2, n_markers);
  values_markers = (double*)mxGetData(outp);
  pargout[0] = outp;

  /* work */
  interpolate_in_triangles(mesh, values, markers, element_id, n_markers, values_markers);
}

