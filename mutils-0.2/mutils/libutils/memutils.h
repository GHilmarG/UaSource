/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef MEMUTILS_H
#define MEMUTILS_H

#include "config.h"
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "debug_defs.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void   inc_memory_usage(size_t size);
  void   dec_memory_usage(size_t size);
  size_t get_thread_memory_usage(void);
  size_t get_total_memory_usage(void);
  void * _mcalloc_threads(const char *var_name, size_t size, size_t start, size_t chunk, int numanode, int thread_id, const char* file, const char* func, const int line);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#define PRINT_MEM_USAGE							\
  printf("\n -- -- MEMORY USAGE %lli K\n\n", memory_usage()/1024)

#ifdef MATLAB_MEX_FILE
#ifdef DEBUG
#pragma message "Using MATLAB memory allocation routines "
#endif
#include "mex.h"
#include "matrix.h"
#define _sys_malloc             malloc
#define _sys_free               free
#define _sys_realloc            realloc
#define _sys_malloc_global      mxMalloc
#define _sys_free_global        mxFree
#define _sys_realloc_global     mxRealloc
#define _sys_persistent(ptr)    mexMakeMemoryPersistent(ptr)
#else
#define _sys_malloc             malloc
#define _sys_free               free
#define _sys_realloc            realloc
#define _sys_malloc_global      malloc
#define _sys_free_global        free
#define _sys_realloc_global     realloc
#define _sys_persistent(ptr)
#endif


#define _mmalloc(func, var, size)					\
  {									\
    void *ptr = 0;							\
    if(size!=0){							\
      ptr = func(size);							\
      if(!ptr){								\
	const char *se = strerror(errno);				\
	EMESSAGE("could not allocate memory (%lli bytes, variable %s)", (long long)(size), #var); \
	EMESSAGE("system error %d: %s", errno, se);			\
	exit(1);							\
      }									\
      inc_memory_usage(size);						\
      if(get_debug_mode()) {						\
	MESSAGE("memory allocated (bytes) %s: %lli (%lli)", #var, (long long)(size), (long long)get_total_memory_usage()); \
	fflush(stdout);							\
      }									\
    } else ptr=NULL;							\
    var = ptr;								\
  }									\


#define _mcalloc(func, var, size)		\
  {						\
    _mmalloc(func, var, size);			\
    if(var) memset(var, 0, size);		\
  }						\


#define _mrealloc(func, var, size, inc)					\
  {									\
    void *ptr;								\
    if(size!=0){							\
      ptr=func(var, size);						\
      if(!ptr){								\
	const char *se = strerror(errno);				\
	EMESSAGE("could not reallocate memory (%lli bytes, variable %s)", (long long)(size), #var); \
	EMESSAGE("system error %d: %s", errno, se);			\
	exit(1);							\
      }									\
      inc_memory_usage(inc);						\
      if(get_debug_mode()) {						\
	MESSAGE("memory reallocated (bytes) %s: %lli (%lli)", #var, (long long)(size), (long long)get_total_memory_usage()); \
	fflush(stdout);							\
      }									\
    } else ptr=NULL;							\
    var = ptr;								\
  }									\


#define _mfree(func, var, size)						\
  {									\
    dec_memory_usage(size);						\
    func(var);								\
    if(get_debug_mode()) {						\
      MESSAGE("memory freed (bytes) %s: %lli (%lli)", #var, (long long)(size), (long long)get_total_memory_usage()); \
      fflush(stdout);							\
    }									\
  }									\



/* 'local' allocation functions, i.e., memory allcated this way  */
/* can not in general be accessed outside of the function that allocates it. */
/* This memory areas can not be used in the MATLAB workspace, but can be used */
/* locally in a MEX function */
#define mmalloc_local(var, size)       _mmalloc(_sys_malloc, var, size);
#define mcalloc_local(var, size)       _mcalloc(_sys_malloc, var, size);
#define mrealloc_local(var, size, inc) _mrealloc(_sys_realloc, var, size, inc);
#define mfree_local(var, size)         _mfree(_sys_free, var, size); var=NULL;


/* 'global' allocation functions, i.e., memory can be used outside */
/* of the function that allocated it provided that it is made persistent */
/* before leaving the function using a call to mpersistent */
#define mmalloc_global(var, size)       _mmalloc(_sys_malloc_global, var, size);
#define mcalloc_global(var, size)       _mcalloc(_sys_malloc_global, var, size);
#define mrealloc_global(var, size, inc) _mrealloc(_sys_realloc_global, var, size, inc);
#define mfree_global(var, size)         _mfree(_sys_free_global, var, size); var=NULL;
#define mpersistent(var)                _sys_persistent(var)


/* Default allocation functions. By default libutils use local memory allocation. */
#define mmalloc(var, size) mmalloc_local(var, size)
#define mcalloc(var, size) mcalloc_local(var, size)
#define mrealloc(var, size, inc) mrealloc_local(var, size, inc)
#define mfree(var, size) mfree_local(var, size)


#define mcalloc_threads(var, size, start, chunk, numanode, thread_id) var = _mcalloc_threads(#var, size, start, chunk, numanode, thread_id, __FILE__, __FUNCTION__, __LINE__);


#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#endif

