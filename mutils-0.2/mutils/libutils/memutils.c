/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "memutils.h"

#ifdef USE_NUMA
#include <numa.h>
#include <numaif.h>
#endif

static size_t total_mem = 0;
static size_t total_thread_mem = 0;
#ifdef USE_OPENMP
#pragma omp threadprivate(total_thread_mem)
#endif /* USE_OPENMP */

void   inc_memory_usage(size_t size)
{
#ifdef USE_OPENMP
#pragma omp atomic
#endif /* USE_OPENMP */
  total_mem += size;
  total_thread_mem += size;
}

void   dec_memory_usage(size_t size)
{
#ifdef USE_OPENMP
#pragma omp atomic
#endif /* USE_OPENMP */
  total_mem -= size;
  total_thread_mem -= size;
}

size_t get_total_memory_usage(void)
{
  return total_mem;
}

size_t get_thread_memory_usage(void)
{
  return total_thread_mem;
}


void *_mcalloc_threads(const char *var_name, size_t size, size_t start, size_t chunk, int numanode, int thread_id, 
		       const char* file, const char* func, const int line)
{
  void *temp = NULL;

#ifdef USE_OPENMP
#pragma omp single copyprivate(temp)
#endif /* USE_OPENMP */
  {
#ifdef WINDOWS
#pragma message ("NUMA binding not implemented on Windows")
    numanode = 0;
#endif
#ifdef _MSC_VER
    temp = malloc(size);
#else
    temp = memalign(4096, size);
#endif
    if(!temp){
      const char *se = strerror(errno);
      EMESSAGE("%s: %s(): %i: could not allocate memory (%lli bytes, variable %s)", file, func, line, (long long)size, var_name);
      EMESSAGE("system error %d: %s", errno, se);
      exit(1);
    }
    inc_memory_usage(size);
  }

#ifdef USE_OPENMP
#pragma omp barrier
#endif /* USE_OPENMP */
  {
    /* MESSAGE("thread %d bind to %x, %x start %d, chunk %d", thread_id, mask, temp, start, chunk); */
#ifdef USE_NUMA
    unsigned long mask = 1 << numanode;
    MESSAGE("mbind %d, %s", mbind((char*)temp+start, chunk, MPOL_BIND, &mask, 3, MPOL_MF_STRICT), strerror(errno));
#endif
    memset((char*)temp+start, 0, chunk);
    /* bzero((char*)temp, size); */
  }

#ifdef USE_OPENMP
#pragma omp barrier
#endif /* USE_OPENMP */

  return temp;
}

/* memutils.c ends here */
