/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef _PARALLEL_H
#define _PARALLEL_H

#include "debug_defs.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef MATLAB_MEX_FILE
#include <matlab/message_id.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  int parallel_set_num_threads(int nthr);
  void parallel_get_info(int *thrid, int *nthr);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _OPENMP_H */
