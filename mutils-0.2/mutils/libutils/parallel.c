/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "parallel.h"
#include <stdlib.h>

int parallel_set_num_threads(int n_omp_threads)
{
#ifdef USE_OPENMP
  char *env_omp_num_threads = NULL;
  env_omp_num_threads = getenv("OMP_NUM_THREADS");
  if(env_omp_num_threads){
    n_omp_threads = atoi(env_omp_num_threads);
  }
  if(!n_omp_threads) n_omp_threads = 1;
  omp_set_num_threads(n_omp_threads);
#else
#ifdef MATLAB_MEX_FILE
  {
    char *env_omp_num_threads = NULL;
    int env_nthr;
    env_omp_num_threads = getenv("OMP_NUM_THREADS");
    if(env_omp_num_threads){
      env_nthr = atoi(env_omp_num_threads);
      if(n_omp_threads>1 || (n_omp_threads==0 && env_nthr>1)){
	USERWARNING("MEX file has been compiled without OpenMP support. Running on 1 CPU", MUTILS_NO_OPENMP);
      }
    }
  }
#endif
  n_omp_threads = 1;
#endif

  return n_omp_threads;
}


void parallel_get_info(int *thrid, int *nthr)
{
#ifdef USE_OPENMP
    *thrid = omp_get_thread_num();
    *nthr  = omp_get_num_threads();
#else
    *thrid = 0;
    *nthr  = 1;
#endif
    DMESSAGE("threads %d/%d", DEBUG_DETAILED, (*thrid), (*nthr));
}
