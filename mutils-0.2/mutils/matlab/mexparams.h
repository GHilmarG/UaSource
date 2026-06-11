/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef _MEXPARAMS_H
#define _MEXPARAMS_H

#include <libutils/config.h>
#include <libutils/mtypes.h>
#include <libutils/debug_defs.h>
#include <mex.h>
#include <matrix.h>

#include "message_id.h"

Double *mex_get_matrix_Double(const mxArray *param, Int *m, Int *n, 
			      const char *varname, const char *sm, const char *sn, int can_be_empty);

mxArray *mex_set_matrix_Double(Double *values, Int m, Int n);

Int *mex_get_matrix_Int(const mxArray *param, Int *m, Int *n, 
			const char *varname, const char *sm, const char *sn, int can_be_empty);

mxArray *mex_set_matrix_Int(Int *values, Int m, Int n);

Int mex_get_integer_scalar(const mxArray *param, const char *varname, int can_be_empty, Int def);


#endif /* _MEXPARAMS_H */
