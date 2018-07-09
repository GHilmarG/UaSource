/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef _MESH_H
#define _MESH_H

#include <libutils/config.h>
#include <libutils/mtypes.h>
#include <libutils/debug_defs.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"

#include "message_id.h"
#endif

/* MESH definition structure */
typedef struct {
  Int          n_elems, n_elem_nodes;
  Int         *elems;
  Int          n_nodes, n_dim;
  Double      *nodes;
  Int          n_neighbors;
  Int         *neighbors;
} t_mesh;


#ifdef MATLAB_MEX_FILE

t_mesh mex2mesh(const mxArray *mesh_struct);

#endif

#endif
