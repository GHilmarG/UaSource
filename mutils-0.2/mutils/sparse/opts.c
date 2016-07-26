/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "opts.h"

#ifdef MATLAB_MEX_FILE

t_opts mex2opts(const mxArray *opts_struct){

  mxArray *field;
  t_opts opts = {0};

  /* default options */
  opts.n_row_entries = -1;
  opts.n_node_dof = 1;
  opts.n_elem_dof = 0;
  opts.symmetric = 0;
  opts.gen_map = 0;

  if(!opts_struct) return opts;

  if(!mxIsStruct(opts_struct)){
    mexErrMsgTxt("opts_struct is not a structure");
  }

  /* n_row_entries */
  field = mxGetField(opts_struct, 0, "n_row_entries");
  opts.n_row_entries = mex_get_integer_scalar(field, "n_row_entries", 1, opts.n_row_entries);

  /* n_node_dof */
  field = mxGetField(opts_struct, 0, "n_node_dof");
  opts.n_node_dof = mex_get_integer_scalar(field, "n_node_dof", 1, opts.n_node_dof);

  /* n_elem_dof */
  field = mxGetField(opts_struct, 0, "n_elem_dof");
  opts.n_elem_dof = mex_get_integer_scalar(field, "n_elem_dof", 1, opts.n_elem_dof);

  /* symmetric */
  field = mxGetField(opts_struct, 0, "symmetric");
  opts.symmetric = mex_get_integer_scalar(field, "symmetric", 1, opts.symmetric);

  /* gen_map */
  field = mxGetField(opts_struct, 0, "gen_map");
  opts.gen_map = mex_get_integer_scalar(field, "gen_map", 1, opts.gen_map);

  return opts;
}

#endif /* MATLAB_MEX_FILE */
