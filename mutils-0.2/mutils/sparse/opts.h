/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef _OPTS_H
#define _OPTS_H

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#include <matrix.h>
#endif

typedef struct {
  int n_row_entries;
  int n_node_dof;
  int n_elem_dof;
  int symmetric;
  int gen_map;
} t_opts;

#ifdef MATLAB_MEX_FILE
#include <matlab/mexparams.h>

t_opts mex2opts(const mxArray *mesh_struct);

#endif

#endif /* _OPTS_H */
