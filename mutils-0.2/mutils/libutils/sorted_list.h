/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef SORTED_LIST_H
#define SORTED_LIST_H

#include "memutils.h"
#include "mtypes.h"

#undef  ALLOC_BLOCKSIZE
#define ALLOC_BLOCKSIZE 16

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void sorted_list_create(Int **list, Int *size);
  void sorted_list_create_double(Double **list, Int *size);
  void sorted_list_create_pair(Int **list, Double **listd, Int *size);

  /* operations on static lists are not inlined to make a more compact code */
  void sorted_list_add(Int **list, Int *nelems, Int *lsize, Int value);
  void sorted_list_add_accum(Int **list, Int *nelems, Int *lsize, Int value, Double **dlist, Double dvalue);


#define SORTED_LIST_LOCATE_H(type)					\
  type sorted_list_locate_##type(type *list, type nelems, type value);

#define SORTED_LIST_LOCATE_C(type)					\
  type sorted_list_locate_##type(type *list,				\
				 type nelems, type value)		\
  {									\
    type l, u;								\
									\
    l = 0;								\
    u = nelems;								\
									\
    /* locate the range by bisection */					\
    /* NOTE: removed since the search is slower for short lists */	\
    /* with which we deal in practice */				\
									\
    /* while(u-l>64){     */						\
    /* m = (l+u)/2;       */						\
    /* if(list[m]>value){ */						\
    /* u=m;		  */						\
    /* } else {		  */						\
    /* l=m;		  */						\
    /* }		  */						\
    /* }		  */						\
									\
    /* locate the value by linear search */				\
    while(l<u){								\
      if(list[l]>=value) break;						\
      l++;								\
    }									\
									\
    return l;								\
  }


#define SORTED_LIST_ADD_STATIC_H(type)				\
  void sorted_list_add_static_##type(type *list, type *nelems,	\
				     type value)

#define SORTED_LIST_ADD_STATIC_C(type)				\
  void sorted_list_add_static_##type(type *list, type *nelems,	\
				     type value)		\
  {								\
    type l, u;							\
								\
    /* locate insert position */				\
    l = sorted_list_locate_##type(list, *nelems, value);	\
    if(l<*nelems && list[l]==value) return;			\
								\
    /* insert into array */					\
    u = *nelems;						\
    while(l!=u){						\
      list[u] = list[u-1];					\
      u--;							\
    }								\
								\
    list[l] = value;						\
    (*nelems)++;						\
  }


#define SORTED_LIST_ADD_STATIC_ACCUM_H(itype, dtype)	\
  void sorted_list_add_static_accum_##itype##_##dtype	\
  (itype *list, itype *nelems,				\
   itype value, dtype *dlist,				\
   dtype dvalue);

#define SORTED_LIST_ADD_STATIC_ACCUM_C(itype, dtype)		\
  void sorted_list_add_static_accum_##itype##_##dtype		\
  (itype *list, itype *nelems,					\
   itype value, dtype *dlist,					\
   dtype dvalue)						\
  {								\
    itype l, u;							\
								\
    /* locate insert position */				\
    l = sorted_list_locate_##itype(list, *nelems, value);	\
    if(l<*nelems && list[l]==value) {				\
      dlist[l] += dvalue;					\
      return;							\
    }								\
								\
    /* insert into array */					\
    u = *nelems;						\
    while(l!=u){						\
      list[u]  = list[u-1];					\
      dlist[u] = dlist[u-1];					\
      u--;							\
    }								\
								\
    list[l]  = value;						\
    dlist[l] = dvalue;						\
    (*nelems)++;						\
  }




#if HAVE_INLINE

  /* define the inline functions */
#ifdef DEBUG
#pragma message "using inline sorted list"
#endif
  static INLINE SORTED_LIST_LOCATE_C(Int)
       static INLINE SORTED_LIST_LOCATE_C(Long)
       static INLINE SORTED_LIST_LOCATE_C(dimType)
       static INLINE SORTED_LIST_LOCATE_C(indexType)

       static INLINE SORTED_LIST_ADD_STATIC_C(Int)
       static INLINE SORTED_LIST_ADD_STATIC_C(dimType)

       static INLINE SORTED_LIST_ADD_STATIC_ACCUM_C(Int, Double)
       static INLINE SORTED_LIST_ADD_STATIC_ACCUM_C(dimType, Double)

#else

  /* declare the functions and define them in the C file */
#ifdef DEBUG
#pragma message "using compiled sorted list"
#endif
  SORTED_LIST_LOCATE_H(Int);
  SORTED_LIST_LOCATE_H(Long);
  SORTED_LIST_LOCATE_H(dimType);
  SORTED_LIST_LOCATE_H(indexType);

  SORTED_LIST_ADD_STATIC_H(Int);
  SORTED_LIST_ADD_STATIC_H(dimType);

  SORTED_LIST_ADD_STATIC_ACCUM_H(Int, Double);
  SORTED_LIST_ADD_STATIC_ACCUM_H(dimType, Double);

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif

