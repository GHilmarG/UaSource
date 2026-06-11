/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef _OSYSTEM_H
#define _OSYSTEM_H

#if defined(_WIN32)
#define WINDOWS
#endif

#ifdef _MSC_VER
/* 'function' : unreferenced inline function has been removed */
#  pragma warning(disable : 4514)
/* 'bytes' bytes padding added after construct 'member_name' */
#  pragma warning(disable : 4820)
/* turn on sse2 on Visual Studio */
#  define __SSE2__
#  ifdef _STDC_
#    define __STRICT_ANSI__
#  endif
#endif /* _MSC_VER */

#ifdef __ICC
/* external function definition with no prior declaration */
#pragma warning disable 1418
/* external declaration in primary source file */
#pragma warning disable 1419
/* floating-point equality and inequality comparisons are unreliable */
#pragma warning disable 1572
/*  #pragma once is obsolete. Use #ifndef guard instead. */
#pragma warning disable 1782
#endif

#ifndef __STRICT_ANSI__
#  define HAVE_INLINE 1
#  ifdef _MSC_VER
#    define INLINE __inline
#  elif defined __GNUC__
#    define INLINE inline
#  else
#    define HAVE_INLINE 0
#ifdef DEBUG
#pragma message "No inline functions - unknown compiler"
#endif
#  endif
#else
#  ifdef __GNUC__
#    define HAVE_INLINE 1
/* in GCC this works despite the ANSI mode */
#    define INLINE __inline__
#  else
#ifdef DEBUG
#pragma message "No inline functions."
#endif
#    define HAVE_INLINE 0
#    define INLINE
#  endif
#endif /* __STRICT_ANSI__ */

#endif /* _OSYSTEM_H */
