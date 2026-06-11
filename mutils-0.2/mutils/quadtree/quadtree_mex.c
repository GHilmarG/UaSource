/*
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

/* libutils headers */
#include <libutils/config.h>
#include <libutils/mtypes.h>
#include <libutils/utils.h>
#include <libutils/sorted_list.h>
#include <libutils/memutils.h>
#include <libutils/tictoc.h>
#include <matlab/mesh.h>
#include <matlab/mexparams.h>
#include <matlab/message_id.h>


/* system headers */
#include <stdlib.h>
#include <stdio.h>
#include <emmintrin.h>
#include <limits.h>
#include <libutils/parallel.h>

#define MAX_TREE_DEPTH (sizeof(Int)*8)
#define ROOT_DEPTH (MAX_TREE_DEPTH-1)
#define MAX_VAL (double)((Int)ULONG_MAX ^ (1<<ROOT_DEPTH))

Double MACHEPS;
#define EMPTY_ELID ((Int)-1)

Int vtk_write2d(char *model_name, Int *elems, Double *nodes, Int *celldata,
		Int nnod, Int nel, Int nnodel);


/* quadrant structure */
#define EMPTY_QUADRANT ((Int)-1)
typedef struct _t_quadrant t_quadrant;
struct  _t_quadrant {
  Uint x_code;
  Uint y_code;
  Uint level;
  size_t parent;
  size_t children[4];
  Int  n_points;        /* how many points are in the quadrant */
  Int  point_id[];      /* point id of the points in quadrant. */
                        /* this thing is actually an array of n_leaf_points point ids */
};


/* quadtree structure - memory area header for the MEX version */
#define QTREE_STR_ID "qTreeMP"
typedef struct {
  char name[8];
  Int n_leaves;
  Int n_quadrants;
  Int n_leaf_points;
  size_t quadrant_size;
  Int n_points;
  Double xmin, xmax;
  Double ymin, ymax;
} t_quadtree;

size_t header_size = sizeof(t_quadtree);

/* memory pool structure */
typedef struct {

  /* complete memory allocated, including the header */
  /* and the subsequent list of quadrants pointed to by base_ptr */
  char *head_ptr;

  /* base_ptr is the pointer to the quadrant array */
  /* However, we can not use t_quadrant since the type definition is incomplete. */
  /* Therefore pointer arithmetic on t_quadrant* is not defined */
  char *base_ptr; 

  /* quadrant_size (size of the t_quadrant structure) */
  /* depends on the n_leaf_points specified at runtime */
  size_t quadrant_size;
  size_t size;
  size_t realloc_size;
  Int ptr;
} t_mempool;


/* global variables */
/* search statistics */
static Long n_leaves             = 1;
static Long n_elems_searched     = 0;
static Double avg_elems_searched = 0;
static Long n_max_elems_searched = 0;


/* lists of elements to be searched */
/* while looking for the element containing a marker */
static Uint nlists   = 0;
static dimType *slist[1024]      = {0};
static size_t  *slist_size[1024] = {0};
static Uint initialized = 0;


/* free list structure */
void quadtree_mex_cleanup(void) {
  Int i;
  
  for(i=0; i<nlists; i++)  {
    if(slist[i]) {
      mfree(slist[i], sizeof(dimType)*slist_size[i][0]);
    }
  }
}


/* TODO fix for larger Int */
INLINE Int pow2roundup (Int x){
    if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}


/* compute the quadrant pointer from the base memory pool address */
/* and quadrant offset */
#define CHILD_POINTER(node, n, mempool)					\
  (t_quadrant*)(mempool->base_ptr + node->children[n])

/* allocate and initialize new leaf quadrant */
/* reallocate memory pool if necessary */
#define CREATE_CHILD(dest, n, _x_code, _y_code, mempool)		\
  {									\
    t_quadrant  *child = NULL;						\
    size_t offset = (char*)dest - mempool->base_ptr;			\
    if(mempool->ptr == mempool->size){					\
      mempool->size += mempool->realloc_size;				\
      mrealloc(mempool->head_ptr,					\
	       header_size + mempool->size*mempool->quadrant_size,	\
	       mempool->realloc_size*mempool->quadrant_size);		\
      mempool->base_ptr = mempool->head_ptr + header_size;		\
      dest = (t_quadrant*)(mempool->base_ptr + offset);			\
    }									\
    dest->children[n] = (size_t)mempool->ptr*mempool->quadrant_size;	\
    mempool->ptr++;							\
    child             = CHILD_POINTER(dest, n, mempool);		\
    child->x_code     = _x_code;					\
    child->y_code     = _y_code;					\
    child->level      = dest->level-1;					\
    child->parent     = offset;						\
    child->children[0]= (size_t)EMPTY_QUADRANT;				\
    child->n_points   = 0;						\
  }									\


t_quadrant *quadtree_locate_codes(t_quadrant *dest, Double nx, Double ny, 
				  t_mempool *mempool)
{
  /* node coordinates out of bounds - point outside of the quadrant */
  Double x = (Double)dest->x_code/MAX_VAL;
  Double y = (Double)dest->y_code/MAX_VAL;
  Double d = 1.0/(Double)(MAX_TREE_DEPTH-dest->level);
  Uint x_code;
  Uint y_code;
  Uint level; 
  Uint bit;
  Uint child;

  if(nx<x || nx>x+d || ny<y || ny>y+d)  return NULL;

  /* fix the case where point is located at the domain boundary */
  if(nx==1.0) nx = nx-MACHEPS;
  if(ny==1.0) ny = ny-MACHEPS;

  x_code = (Uint)(nx*MAX_VAL);
  y_code = (Uint)(ny*MAX_VAL);
  level = dest->level-1;

  bit = 1 << level;
  while((dest)->children[0] != EMPTY_QUADRANT){
    child  = ((x_code & bit)>>level); level--;
    child += ((y_code & bit)>>level);
    dest = CHILD_POINTER(dest, child, mempool);
    bit >>= 1;
  }

  return dest;
}


/* Incrementally build a quadtree from nodes. */
/* Add nodes in sequence, quadtree structure is refined in the process. */
/* Internally the wuadtree structure is built from a normalized domain, */
/* i.e., coordinates from [0,1]. The coordinates of added points */\
/* are normalized as we go. */
void quadtree_add_node(t_quadrant *dest, Double *nodes, Int point_id, Int n_leaf_points,
		       Int *n_qtree_points, t_mempool *mempool)
{
  Double nx = nodes[point_id*2+0];
  Double ny = nodes[point_id*2+1];
  
  {
    /* normalize coordinates */
    t_quadtree *tree = (t_quadtree *)mempool->head_ptr;
    nx = (nx - tree->xmin)/(tree->xmax - tree->xmin);
    ny = (ny - tree->ymin)/(tree->ymax - tree->ymin);
  }

  /* look for the destination quadrant */
  if(dest->children[0] != EMPTY_QUADRANT){
    dest = quadtree_locate_codes(dest, nx, ny, mempool);
  }

  /* nothing to do - point outside of the quadrant domain */
  if(!dest) return;

  /* safequard - quadtree maximum level exceeded */
  if(dest->level==0) return;

  /* if the quadrant has free space simply add the node */
  if(dest->n_points < n_leaf_points){
    dest->point_id[dest->n_points++] = point_id;
    (*n_qtree_points)++;
    return;
  }

  /* split the leaf (dest) and reassign the nodes to new quadtree leaves */

  /* do not clear the node information in the parent */
  /* useful when leaves are empty and we want to have */
  /* some information about nearby points/elements during search */

  {
    Uint bit = 1<<(dest->level-1);
    CREATE_CHILD(dest, 0, (dest->x_code)      , (dest->y_code)      , mempool);
    CREATE_CHILD(dest, 1, (dest->x_code) | bit, (dest->y_code)      , mempool);
    CREATE_CHILD(dest, 2, (dest->x_code)      , (dest->y_code) | bit, mempool);
    CREATE_CHILD(dest, 3, (dest->x_code) | bit, (dest->y_code) | bit, mempool);
  }

  /* update number of leaves */
  n_leaves += 3;
  
  /* add the old nodes directly to the correct child quadrant */
  {
    Int ptid;

    t_quadrant *child;
    Uint level = dest->level-1;
    Uint bit   = 1 << level;
    Uint x_code, y_code;

    for(ptid=0; ptid<n_leaf_points; ptid++){

      nx = nodes[dest->point_id[ptid]*2+0];
      ny = nodes[dest->point_id[ptid]*2+1];

      {
	/* normalize coordinates */
	/* NOTE: memory pool might have been reallocated, refresh tree pointer! */
	t_quadtree *tree = (t_quadtree *)mempool->head_ptr;
	nx = (nx - tree->xmin)/(tree->xmax - tree->xmin);
	ny = (ny - tree->ymin)/(tree->ymax - tree->ymin);
      }

      x_code = (Uint)(nx*MAX_VAL);
      y_code = (Uint)(ny*MAX_VAL);

      child = CHILD_POINTER(dest, 
			    ((x_code & bit)>>level) + 
			    ((y_code & bit)>>(level-1)), mempool);

      child->point_id[child->n_points++] = dest->point_id[ptid];
    }
  }

  /* add the new node recursively */
  quadtree_add_node(dest, nodes, point_id, n_leaf_points, n_qtree_points, mempool);
}


/* linearize the quadtree - extract leaves in Z-curve ordering */
#ifdef _MSC_VER
#pragma auto_inline(off)
#endif
void quadtree_extract_leaves(t_quadrant *dest, t_quadrant **tree_leaves, Int *itree_leaves, t_mempool *mempool)
{

  /* store the leaves */
  if(dest->children[0]==EMPTY_QUADRANT){
    tree_leaves[*itree_leaves] = dest;
    (*itree_leaves)++;
    return;
  }

  /* traverse the subtrees */
  quadtree_extract_leaves(CHILD_POINTER(dest, 0, mempool), tree_leaves, itree_leaves, mempool);
  quadtree_extract_leaves(CHILD_POINTER(dest, 1, mempool), tree_leaves, itree_leaves, mempool);
  quadtree_extract_leaves(CHILD_POINTER(dest, 2, mempool), tree_leaves, itree_leaves, mempool);
  quadtree_extract_leaves(CHILD_POINTER(dest, 3, mempool), tree_leaves, itree_leaves, mempool);
}


/* linearize the quadtree - extract points in Z-curve ordering */
#ifdef _MSC_VER
#pragma auto_inline(off)
#endif
void quadtree_extract_points(t_quadrant *dest, Int *points, Int *points_ptr, t_mempool *mempool)
{
  Int i;

  /* copy point data from the leaves */
  if(dest->children[0]==EMPTY_QUADRANT){
    if(dest->n_points){
      for(i=0; i<dest->n_points; i++){
	points[(*points_ptr)+i] = dest->point_id[i]+1;
      }
      (*points_ptr) += dest->n_points;
    }
    return;
  }

  /* traverse the subtrees */
  quadtree_extract_points(CHILD_POINTER(dest, 0, mempool), points, points_ptr, mempool);
  quadtree_extract_points(CHILD_POINTER(dest, 1, mempool), points, points_ptr, mempool);
  quadtree_extract_points(CHILD_POINTER(dest, 2, mempool), points, points_ptr, mempool);
  quadtree_extract_points(CHILD_POINTER(dest, 3, mempool), points, points_ptr, mempool);
}


#ifndef __SSE2__
#pragma message ("Using non-SSE vector product.")

#ifdef __STRICT_ANSI__
#define cross(p, a, b)					\
  ((b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]))

#else
INLINE double cross(const double *p, const double *a, const double *b)
{
  return ((b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]));
}
#endif /* __STRICT_ANSI__ */

#else

#ifdef DEBUG
#pragma message ("Using SSE vector product.")
#endif

/* TODO: single precision */
INLINE double cross(const double *p, const double *a, const double *b)
{
  __m128d v_p , v_a , v_b;
  double temp;

  /* load data to registers */
  v_p = _mm_load_pd(p);
  v_a = _mm_load_pd(a);
  v_b = _mm_load_pd(b);

  v_b = _mm_sub_pd(v_b, v_a); /* b-a */
  v_p = _mm_sub_pd(v_p, v_a); /* p-a */

  v_p = _mm_mul_pd(v_b, _mm_shuffle_pd(v_p, v_p, _MM_SHUFFLE2 (0,1)));

  _mm_store_sd (&temp, _mm_sub_pd(v_p, _mm_shuffle_pd(v_p, v_p, _MM_SHUFFLE2 (0,1))));
  return temp; 
}
#endif /* __SSE2__ */



#define ENQUEUE_NEIGHBOR(n)						\
  if((n)!=EMPTY_ELID){							\
    if(thr_slist_nel==thr_slist_size){					\
      mrealloc(slist[thrid],sizeof(dimType)*2*thr_slist_size,sizeof(dimType)*thr_slist_size); \
      slist_size[thrid][0] *= 2;					\
      thr_slist = slist[thrid];						\
      thr_slist_size = slist_size[thrid][0];				\
    }									\
    thr_slist[thr_slist_nel++] = (n);					\
  }									\


int quadtree_locate_element(Int elid, Int marker_id, Double *markerX, t_mesh mesh, 
			    Int *map, Long *nel_searched, Int thrid)
{
  dimType n1, n2, n3;
  const Double *a, *b, *c;

  /* lists of elements to be searched */
  /* while looking for the element containing a marker */
  dimType *thr_slist     = slist[thrid];
  size_t thr_slist_size = slist_size[thrid][0];
  size_t thr_slist_nel  = 0;
  size_t thr_slist_ptr  = 0;
  size_t li;

  *nel_searched = 0;

  /* search elements and their neighbors */
  while(1){

    /* has the point in triangle test been performed? */
    if((elid != EMPTY_ELID) && (map[elid] != marker_id)) {

      /* perform the point in triangle test */
      map[elid] = marker_id;
      (*nel_searched)++;


      li = (size_t)elid*mesh.n_elem_nodes;
      a = mesh.nodes + (size_t)2*(mesh.elems[li+0]-ONE_BASED_INDEX);
      b = mesh.nodes + (size_t)2*(mesh.elems[li+1]-ONE_BASED_INDEX);
      c = mesh.nodes + (size_t)2*(mesh.elems[li+2]-ONE_BASED_INDEX);


      li = (size_t)elid*mesh.n_neighbors;
      n1 = mesh.neighbors[li + 2]-ONE_BASED_INDEX;
      n2 = mesh.neighbors[li + 0]-ONE_BASED_INDEX;
      n3 = mesh.neighbors[li + 1]-ONE_BASED_INDEX;

      /* Point in triangle, half-planes test. */
      /* Relies on counter-clockwise ordering of triangle nodes */
      /* and on a correct order of triangle neighbors, i.e., */
      /* neighbor 1 lies accross the edge opposite to node 1, and so on. */

      /* Searching for the containing triangle is done using */
      /* Green and Sibson algorithm. Termination of the algorithm is assured */
      /* through a test map (every triangle is only tested once) */
      /* and a queue of triangles to be checked if the simple approach fails. */
      /* In the worst-case scenario all elements are verified. */

      if(cross(markerX,a,b) < 0){
	elid = n1;
	ENQUEUE_NEIGHBOR(n2);
	ENQUEUE_NEIGHBOR(n3);
	continue;
      }

      if(cross(markerX,b,c) < 0){
	elid = n2;
	ENQUEUE_NEIGHBOR(n1);
	ENQUEUE_NEIGHBOR(n3);
	continue;
      }

      if(cross(markerX,c,a) < 0){
	elid = n3;
	ENQUEUE_NEIGHBOR(n1);
	ENQUEUE_NEIGHBOR(n2);
	continue;
      }

      return elid;
    }

    /* anything left in the queue? */
    do{
      if(thr_slist_ptr == thr_slist_nel) return -1;
      elid = thr_slist[thr_slist_ptr++];
    } while(elid==EMPTY_ELID);
  }

  /* TODO: if no element is found the procedure takes too long time: */
  /* all elements are checked! */
}


/***********************************************************/
/*              MATLAB INTERFACE                           */
/***********************************************************/

mxArray *qtree2mex(t_quadtree *tree, size_t tree_size){
#define n_fieldnames 5
  const char *fieldnames[n_fieldnames] = {"qtree", "n_leaves", "n_leaf_points", "n_quadrants", "n_points"};
  mxArray *outp = mxCreateStructMatrix(1, 1, n_fieldnames, fieldnames);
  mxArray *field;
  Int n = 0;
  
  /* TODO: fix for larger Int size */
  /* field = set_pointer_handle(tree, "QuadTree"); */
  field = mxCreateNumericMatrix(0, 0, mxUINT8_CLASS,mxREAL);
  mxSetData(field, (void*)tree);
  mxSetN(field, 1);
  mxSetM(field, tree_size);
  mxSetField(outp, 0, fieldnames[n++], field);
  
  field = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  ((Int*)mxGetData(field))[0] = tree->n_leaves;
  mxSetField(outp, 0, fieldnames[n++], field);
  
  field = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  ((Int*)mxGetData(field))[0] = tree->n_leaf_points;
  mxSetField(outp, 0, fieldnames[n++], field);
  
  field = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  ((Int*)mxGetData(field))[0] = tree->n_quadrants;
  mxSetField(outp, 0, fieldnames[n++], field);

  field = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  ((Int*)mxGetData(field))[0] = tree->n_points;
  mxSetField(outp, 0, fieldnames[n++], field);

  return outp;
}


t_quadtree* mex2qtree(const mxArray *qtree_struct){

  mxArray *field;
  t_quadtree *qtree;

  if(!mxIsStruct(qtree_struct)){
    USERERROR("qtree_struct is not a structure", MUTILS_INVALID_PARAMETER);
  }

  /* quadtree memory pointer */
  field = mxGetField(qtree_struct, 0, "qtree");
  if(!field){
    USERERROR("qtree_struct is not a valid quadtree", MUTILS_INVALID_PARAMETER);
  }

  qtree = (t_quadtree*)mxGetData(field);

  /* verify the contents - memory area header */
  qtree->name[7] = 0;

  if(strcmp(qtree->name, QTREE_STR_ID)){
    USERERROR("qtree_struct is not a valid quadtree - invalid header", MUTILS_INVALID_PARAMETER);
  }

  return qtree;
}


mxArray *mex_quadtree_create(int nargin, const mxArray *pargin[])
{
  Double *points = NULL;
  Int n_points;
  Int n_dim;
  Int arg = 1, i;
  Int n_leaf_points;

  /* domain size */
  Double   xmin=0, xmax=1, ymin=0, ymax=1;
  t_mempool mempool = {0};
  size_t initial_size = 0;
  t_quadtree *qtree = NULL;
  t_quadrant *root  = NULL;
  Int n_qtree_points = 0;

  if(!initialized){
    initialized = 1;
    mexAtExit(quadtree_mex_cleanup);
  }

  if(nargin<2){
    USERERROR("Usage: quadtree = quadtree('create', POINTS, [xmin], [xmax], [ymin], [ymax], [max_points_in_quadrant])", MUTILS_INVALID_PARAMETER);
  }

  /* POINTS */
  n_dim = 2;
  n_points = 0;
  points = mex_get_matrix_Double(pargin[arg++], &n_dim, &n_points, "POINTS", "2", "number of points", 0);

  /* domain extents */
  i = 1;
  if(nargin>arg){
    xmin = mex_get_matrix_Double(pargin[arg++], &i, &i, "xmin", "1", "1", 0)[0];
  } 
    
  if(nargin>arg){
    xmax = mex_get_matrix_Double(pargin[arg++], &i, &i, "xmax", "1", "1", 0)[0];
  } 
    
  if(nargin>arg){
    ymin = mex_get_matrix_Double(pargin[arg++], &i, &i, "ymin", "1", "1", 0)[0];
  } 
    
  if(nargin>arg){
    ymax = mex_get_matrix_Double(pargin[arg++], &i, &i, "ymax", "1", "1", 0)[0];
  } 
    
  /* maximum number of points in quadrant */
  if(nargin>arg){
    n_leaf_points = mex_get_integer_scalar(pargin[arg++], "max_points_in_quadrant", 0, 0);
    arg++;
  } else {
    n_leaf_points = 1;
  }
  if(n_leaf_points>n_points) n_leaf_points = n_points;
  if(n_leaf_points<1) n_leaf_points = 1;


  /* setup the memory pool */
  /* Allocate roughly the correct amount of memory */
  /* for the case when points are spread uniformly in space. */
  initial_size = (size_t)n_points*2/pow2roundup(n_leaf_points);
  mempool.size = initial_size;
  mempool.realloc_size = mempool.size/2;
  mempool.ptr  = 1;
  mempool.quadrant_size = sizeof(t_quadrant) + sizeof(Int)*n_leaf_points;
  mcalloc(mempool.head_ptr, header_size + mempool.size*mempool.quadrant_size);
  mempool.base_ptr = mempool.head_ptr + header_size;

  /* set real domain dimensions for coordinate normalization */
  qtree = (t_quadtree*)mempool.head_ptr;
  qtree->xmin = xmin;
  qtree->xmax = xmax;
  qtree->ymin = ymin;
  qtree->ymax = ymax;

  /* add root quadrant */
  root = (t_quadrant*)mempool.base_ptr;
  root->level  = ROOT_DEPTH;
  root->x_code = 0;
  root->y_code = 0;
  root->n_points = 0;
  root->children[0] = (size_t)EMPTY_QUADRANT;
  root->parent = (size_t)EMPTY_QUADRANT;
  n_leaves = 1;

  /* run */
  n_qtree_points = 0;
  for(i=0; i<n_points; i++){

    /* memory pool can be reallocated in quadtree_add_node */
    root = (t_quadrant*)mempool.base_ptr;
    quadtree_add_node(root, points, i, n_leaf_points, &n_qtree_points, &mempool);
  }

  /* fill the memory header */
  qtree = (t_quadtree*)mempool.head_ptr;
  strncpy(qtree->name, QTREE_STR_ID, 8);
  qtree->n_leaves       = n_leaves;
  qtree->n_quadrants    = mempool.ptr;
  qtree->n_leaf_points  = n_leaf_points;
  qtree->quadrant_size  = mempool.quadrant_size;
  qtree->n_points       = n_qtree_points;
  
  if(n_qtree_points != n_points){
    MESSAGE("Some of the points were outside of the specified domain:\n\n\
    (xmin=%.1e, xmax=%.1e, ymin=%.1e, ymax=%.1e)\n\nand were not added to the quadtree.\n \
Please specify correct domain extents.", xmin, xmax, ymin, ymax);
  }

  /* reallocate memory using MATLAB's allocation routines */
  {
    t_quadtree *_qtree;
    mmalloc_global(_qtree, header_size + mempool.size*mempool.quadrant_size);
    memcpy(_qtree, qtree, header_size + mempool.size*mempool.quadrant_size);
    free(qtree);
    qtree = _qtree;
    mpersistent(qtree);
  }

  return qtree2mex(qtree, header_size + mempool.size*mempool.quadrant_size);
}


mxArray *mex_quadtree_locate(int nargin, const mxArray *pargin[])
{
  Int arg = 1;
  Int *element_map = NULL;
  Int n_dim, temp;
  t_quadtree *tree = NULL;
  t_mempool mempool = {0};
  t_mesh mesh = {0};
  Int      n_markers;
  Double  *markers;
  Int     *elids = NULL; 
  mxArray *outp = NULL;
  Int inplace = 0;
   
  if(!initialized){
    initialized = 1;
    mexAtExit(quadtree_mex_cleanup);
  }
    
  if(nargin<4){
    USERERROR("Usage: [MAP, stats] = quadtree('locate', quadtree, MESH, MARKERS, [MAP], [inplace])", MUTILS_INVALID_PARAMETER);
  }
  
  /* qtree structure */
  tree                   = mex2qtree(pargin[arg++]);
  mempool.head_ptr       = (char*)tree;
  mempool.base_ptr       = mempool.head_ptr + header_size;
  mempool.quadrant_size  = tree->quadrant_size;
  mempool.ptr            = tree->n_quadrants;
  
  /* triangular mesh structure */
  mesh = mex2mesh(pargin[arg++]);
  if(!mesh.neighbors){
    USERERROR("MESH must contain NEIGHBORS information", MUTILS_INVALID_MESH);
    return NULL;
  }

  /* MARKERS */
  n_dim = 2;
  n_markers = 0;
  markers = mex_get_matrix_Double(pargin[arg++], &n_dim, &n_markers, "MARKERS", "2", "number of markers", 0);

  /* optional - previous marker-to-element map to use */
  if(nargin>=5){
    temp = 1;
    element_map = mex_get_matrix_Int(pargin[arg++], &temp, &n_markers, "MAP", "1", "number of markers", 1);
  }

  /* optional - inplace map. Existing MAP input will be overwritten and returned as output. */
  /* Not allowed in MATLAB, so be careful and make sure MAP is not used elsewhere/linked to. */
  if(nargin>=6 && element_map){
    inplace = mex_get_integer_scalar(pargin[arg++], "inplace", 1, 0);
  }

  if(inplace){
    outp  = (mxArray *)pargin[4];
    elids = element_map;
  }

  /* MEX output, needs to be global and persistent */
  if(!outp){
    mcalloc_global(elids, sizeof(Int)*n_markers);
    mpersistent(elids);
  }
  n_elems_searched     = 0;
  n_max_elems_searched = 0;
  
  /* use default/environment defined number of threads */
  parallel_set_num_threads(0);

#ifdef USE_OPENMP
#pragma omp parallel
#endif
  {
    Ulong i;
    int thrid, nthr;
    Int marker_start, marker_end;
    Int elid = EMPTY_QUADRANT;
    Long nel_searched;
    t_quadrant *quadrant = NULL;
    Int thr_n_elems_searched = 0;
    Int thr_n_max_elems_searched = 0;
    t_quadrant *root  = NULL;
    Int blk_size;

    Int *map;
    mcalloc(map,   sizeof(Int)*mesh.n_elems);
  
    /* locate the markers in the elements using the quadtree */
    root     = (t_quadrant*)mempool.base_ptr;
    
    parallel_get_info(&thrid, &nthr);

    blk_size     = n_markers/nthr+1;
    marker_start = thrid*blk_size;
    marker_end   = (thrid+1)*blk_size;
    marker_end   = MIN(n_markers, marker_end);

    /* global list initialization */
    nlists = MAX(nlists, nthr);
    if(slist[thrid]==NULL){

      /* allocate a lot to avoid page sharing between threads */
      mmalloc(slist[thrid], sizeof(dimType)*4096);
      mmalloc(slist_size[thrid], sizeof(size_t)*4096);
      slist_size[thrid][0] = 4096;
    }

    for(i=marker_start; i<marker_end; i++){
    
      /* prefetch markers - non-temporal to make space */
      /* for the qtree structure in the CPU caches */
      /* if(i+16<marker_end) _mm_prefetch(((char*)markers)+(i+16), _MM_HINT_NTA); */

      if(element_map){
	/* if(i+16<marker_end) _mm_prefetch(((char*)element_map)+(i+16), _MM_HINT_NTA); */
      	elid = element_map[i] - ONE_BASED_INDEX;
      }

      if(!element_map || elid==EMPTY_QUADRANT){

      	/* Locate the quadrant. */
      	/* quadrant is needed only to get some 'nearby' element id. */
      	/* The correct element containing the marker is located */
      	/* by searching the element neighbors. */

      	/* uptree traversal does not speed up things at all */
      	/* even if the input points are reasonably sorted   */
      	/* quadrant = quadtree_locate_sorted(quadrant, markers[i*2+0], markers[i*2+1], &mempool); */

      	Double nx, ny;

      	/* normalize coordinates */
      	nx = (markers[i*2+0] - tree->xmin)/(tree->xmax - tree->xmin);
      	ny = (markers[i*2+1] - tree->ymin)/(tree->ymax - tree->ymin);
      	quadrant = quadtree_locate_codes(root, nx, ny, &mempool);

      	elid = EMPTY_QUADRANT;

      	if(quadrant){

      	  /* Find a nearby node in the quadtree. */
      	  /* If the given quadrant is empty, */
      	  /* return a node stored in first non-empty parent */
      	  while(1){
	  
      	    /* no elements in this quadrant - go up the tree */
      	    if(quadrant->n_points == 0){
      	      if(quadrant->parent == EMPTY_QUADRANT) break;
      	      quadrant = (t_quadrant*)(quadrant->parent + mempool.base_ptr);
      	      continue;
      	    } else {
      	      elid = quadrant->point_id[0];
      	      break;
      	    }
      	  }
      	}

      	if(elid == EMPTY_QUADRANT){

      	  /* not found. */
      	  elids[i] = ONE_BASED_INDEX + elid;
      	  continue;
      	}
      }

      /* find containing element */
      /* NOTE: coordinate normalization not needed here since we do a mesh search, */
      /* not a quadtree search. */
      elid     = quadtree_locate_element(elid, i+1, markers+i*2, mesh, map, &nel_searched, thrid);
      elids[i] = ONE_BASED_INDEX + elid;
      thr_n_elems_searched += nel_searched;
      thr_n_max_elems_searched = MAX(thr_n_max_elems_searched, nel_searched);
    }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
    n_elems_searched += thr_n_elems_searched;

#ifdef USE_OPENMP
#pragma omp critical
#endif
    n_max_elems_searched = MAX(n_max_elems_searched, thr_n_max_elems_searched);

    mfree(map,   sizeof(Int)*mesh.n_elems);
  }
    
  avg_elems_searched = (Double)n_elems_searched/n_markers;

  if(!outp) outp = mex_set_matrix_Int(elids, 1, n_markers);

  return outp;
}


mxArray *mex_quadtree_reorder(int nargin, const mxArray *pargin[])
{
  Int arg = 1;
  t_quadtree *tree = NULL;
  t_mempool mempool = {0};
  Int *I;
  Int points_ptr = 0;

  mxArray *outp;

  if(!initialized){
    initialized = 1;
    mexAtExit(quadtree_mex_cleanup);
  }
  
  if(nargin<2){
    USERERROR("Usage: I = quadtree('reorder', quadtree)", MUTILS_INVALID_PARAMETER);
  }

  tree                   = mex2qtree(pargin[arg++]);
  mempool.head_ptr       = (char*)tree;
  mempool.base_ptr       = mempool.head_ptr + header_size;
  mempool.quadrant_size  = tree->quadrant_size;
  mempool.ptr            = tree->n_quadrants;

  /* MEX output, needs to be global and persistent */
  mmalloc_global(I, sizeof(Int)*tree->n_points);
  mpersistent(I);

  /* extract nodes in the Z-ordering */
  quadtree_extract_points((t_quadrant*)mempool.base_ptr, I, &points_ptr, &mempool);
  
  outp = mex_set_matrix_Int(I, 1, points_ptr);

  return outp;
}


void mex_vtkwrite(int nargin, const mxArray *pargin[])
{
  Ulong i;
  
  /* prepare vtk data */
  t_quadrant **tree_leaves;
  Int itree_leaves = 0;
  Int arg = 1;

  t_quadtree *tree;
  t_mempool mempool      = {0};
  t_quadrant *root;

  double *vtk_nodes;
  Int    *vtk_elems;
  Int    *vtk_celld;
  Int     n_cells = 0;
  char    fname[512];

  if(!initialized){
    initialized = 1;
    mexAtExit(quadtree_mex_cleanup);
  }

  if(nargin<2){
    USERERROR("Usage: quadtree('vtkwrite', quadtree, [file_name])", MUTILS_INVALID_PARAMETER);
  }

  /* quadtree */
  tree                   = mex2qtree(pargin[arg++]);
  mempool.head_ptr       = (char*)tree;
  mempool.base_ptr       = mempool.head_ptr + header_size;
  mempool.quadrant_size  = tree->quadrant_size;
  mempool.ptr            = tree->n_quadrants;

  /* file name */
  if(nargin>2){
    if(!mxIsChar(pargin[arg])) USERERROR("'file_name' must be a string", MUTILS_INVALID_PARAMETER);
    if(0!=mxGetString(pargin[arg], fname, 511)) USERERROR("file_name too long, can be maximum 511 characters.", MUTILS_INVALID_PARAMETER);
  } else {
    sprintf(fname, "%s", "quadtree");
  }

  root  = (t_quadrant*)mempool.base_ptr;
  mcalloc(tree_leaves, sizeof(t_quadrant*)*n_leaves);
  quadtree_extract_leaves(root, tree_leaves, &itree_leaves, &mempool);

  mcalloc(vtk_nodes, sizeof(double)*4*n_leaves*2);
  mcalloc(vtk_elems, sizeof(Int)*4*n_leaves);
  mcalloc(vtk_celld, sizeof(Int)*1*n_leaves);

  for(i=0; i<n_leaves; i++){

    double mix, miy, max, may;
    double dx, dy;

    dx  = (1<<(tree_leaves[i]->level))/MAX_VAL;
    dy  = (1<<(tree_leaves[i]->level))/MAX_VAL;
    mix = tree_leaves[i]->x_code/MAX_VAL;
    miy = tree_leaves[i]->y_code/MAX_VAL;
    max = dx + mix;
    may = dy + miy;
    dx  = dx*0.05;
    dy  = dy*0.05;
    
    vtk_nodes[i*4*2 + 0*2 + 0] = mix+dx;
    vtk_nodes[i*4*2 + 0*2 + 1] = miy+dy;

    vtk_nodes[i*4*2 + 1*2 + 0] = max-dx;
    vtk_nodes[i*4*2 + 1*2 + 1] = miy+dy;

    vtk_nodes[i*4*2 + 2*2 + 0] = max-dx;
    vtk_nodes[i*4*2 + 2*2 + 1] = may-dy;

    vtk_nodes[i*4*2 + 3*2 + 0] = mix+dx;
    vtk_nodes[i*4*2 + 3*2 + 1] = may-dy;

    vtk_elems[i*4 + 0] = i*4 + 0;
    vtk_elems[i*4 + 1] = i*4 + 1;
    vtk_elems[i*4 + 2] = i*4 + 2;
    vtk_elems[i*4 + 3] = i*4 + 3;

    vtk_celld[i] = tree_leaves[i]->n_points; 
    n_cells += vtk_celld[i];
  }

  vtk_write2d(fname, vtk_elems, vtk_nodes, vtk_celld, n_leaves*4, n_leaves, 4);
  
  mfree(tree_leaves, sizeof(t_quadrant*)*n_leaves);
  mfree(vtk_nodes, sizeof(double)*4*n_leaves*2);
  mfree(vtk_elems, sizeof(int)*4*n_leaves);
  mfree(vtk_celld, sizeof(int)*1*n_leaves);
}


Int vtk_write2d(char *model_name, Int *elems, Double *nodes, Int *celldata,
		Int nnod, Int nel, Int nnodel)
{
  FILE *out_vtk;
  long i;
  char file_name[512+4];

  sprintf(file_name, "%s.vtk", model_name);
  out_vtk=fopen(file_name, "w");

  fprintf(out_vtk,"# vtk DataFile Version 3.0\n");
  fprintf(out_vtk,"my cool data\n");
  fprintf(out_vtk,"ASCII\n");
  fprintf(out_vtk,"DATASET UNSTRUCTURED_GRID\n");

  fprintf(out_vtk,"POINTS %d double\n", nnod);
  for (i=0;i<nnod;i++){
    fprintf(out_vtk,"%lf %lf 0.0\n", nodes[2*i+0], nodes[2*i+1]);
  }

  fprintf(out_vtk,"CELLS %d %d\n", nel, (1+4)*nel);
  for (i=0;i<nel;i++){
    fprintf(out_vtk,"4 %d %d %d %d\n",
	    elems[nnodel*i+0], elems[nnodel*i+1], elems[nnodel*i+2], elems[nnodel*i+3]);
  }
  fprintf(out_vtk,"CELL_TYPES %d\n", nel);
  for (i=0;i<nel;i++){
    fprintf(out_vtk,"9\n");
  }


  fprintf(out_vtk,"CELL_DATA %d\n", nel);
  fprintf(out_vtk,"SCALARS n_nodes_in_quadrant long 1\n");
  fprintf(out_vtk,"LOOKUP_TABLE default\n");
  for (i=0;i<nel;i++){
    fprintf(out_vtk,"%d\n", (int)celldata[i]);  
  }

  fclose(out_vtk); 
  return 0;
}


mxArray *mex_quadtree_stats(void)
{
#undef n_fieldnames
#define n_fieldnames 4
  const char *fieldnames[n_fieldnames] = {"n_elems_searched", "avg_elems_searched", "n_max_elems_searched", "list_size"};
  mxArray *outp = mxCreateStructMatrix(1, 1, n_fieldnames, fieldnames);
  mxArray *field;
  Int n = 0;

  /* TODO fix for larger Int and smaller Double */
  field = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
  ((Long*)mxGetData(field))[0] = n_elems_searched;
  mxSetField(outp, 0, fieldnames[n++], field);
  
  field = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  ((Double*)mxGetData(field))[0] = avg_elems_searched;
  mxSetField(outp, 0, fieldnames[n++], field);

  field = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
  ((Long*)mxGetData(field))[0] = n_max_elems_searched;
  mxSetField(outp, 0, fieldnames[n++], field);

  field = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  if(nlists){
    ((Int*)mxGetData(field))[0] = slist_size[0][0];
  } else {
    ((Int*)mxGetData(field))[0] = -1;
  }
  mxSetField(outp, 0, fieldnames[n++], field);

  return outp;
}


void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[])
{
  int arg = 0;
  char cmd[256];

  /* get machine epsilon */
  MACHEPS = macheps();
 
  if (nargin < 1){
    USERERROR("Wrong usage. Type 'help quadtree' for examples.", MUTILS_INVALID_PARAMETER);
  }

  /* command */
  {
    if(!mxIsChar(pargin[arg])){
      USERERROR("command parameter must be a string", MUTILS_INVALID_PARAMETER);
    }
    mxGetString(pargin[arg], cmd, 255);
    arg++;
  }

  if(!strcmp(cmd, "create")){
    if(nargout>0){
      pargout[0] = mex_quadtree_create(nargin, pargin);
    }
    return;
  }

  if(!strcmp(cmd, "vtkwrite")){
    mex_vtkwrite(nargin, pargin);
    return;
  }

  if(!strcmp(cmd, "locate")){
    if(nargout>0){
      pargout[0] = mex_quadtree_locate(nargin, pargin);
    }
    if(nargout>1){
      pargout[1] = mex_quadtree_stats();
    }
    return;
  }

  if(!strcmp(cmd, "reorder")){
    if(nargout>0){
      pargout[0] = mex_quadtree_reorder(nargin, pargin);
    }
    return;
  }

  USERERROR("unknown command", MUTILS_INVALID_PARAMETER);
}
