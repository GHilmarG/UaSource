/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef MTYPES_H
#define MTYPES_H

#include "config.h"
#include <limits.h>


/* #ifdef linux */
/* #include <bits/types.h> */
/* #else */
/* #include <sys/types.h> */
/* typedef uint64_t __uint64_t; */
/* typedef uint32_t __uint32_t; */
/* #endif */

/* typedef double mfloat; */
/* typedef double mdouble; */

/* typedef int              mint; */
/* typedef unsigned         muint; */
/* typedef unsigned long    muintl; */

/* #if defined (__x86_64__) || defined (__64BIT__) */
/* typedef __uint64_t       mIndexType; */
/* typedef __uint32_t       mDimType; */
/* #define MAX_DIM_TYPE     0xffffffff */
/* #else */
/* typedef __uint32_t       mIndexType; */
/* typedef __uint32_t       mDimType; */
/* #define MAX_DIM_TYPE     0xffffffff */
/* #endif */



typedef double tfloat;
typedef unsigned long tu64;
typedef long ti64;
typedef unsigned int tu32;
typedef int ti32;

typedef int Int;
#define MaxInt UINT_MAX
typedef unsigned int Uint;
#define MaxUint UINT_MAX
typedef long Long;
#define MaxLong LONG_MAX
typedef unsigned long Ulong;
#define MaxUlong ULONG_MAX
typedef double Double;
typedef unsigned long MemoryOffset;

#define IS_TYPE_SIGNED(type)  (((type)-1)<0)

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#define ONE_BASED_INDEX 1
typedef mwIndex       indexType;
typedef mwSize          dimType;
#else
#define ONE_BASED_INDEX 0
typedef unsigned long indexType;
typedef unsigned        dimType;
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  tfloat macheps(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#define managed_type_cast(type, var, val, errmsg)			\
  {									\
    var = (type)val;							\
    if(var!=val){							\
      USERERROR("%s", MUTILS_INCOMPATIBLE_TYPES, errmsg);		\
    }									\
  }
  
#endif /*  _MTYPES_H */
