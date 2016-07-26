/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "utils.h"
#include "mtypes.h"

#undef  ALLOC_BLOCKSIZE
#define ALLOC_BLOCKSIZE 16

#define CREATE_SORTED_LIST_H(type)					\
  typedef struct {							\
    type *ptr;								\
    Int  size;								\
    Int  nelems;							\
  } t_list_##type;							\
									\
  t_list_##type sorted_list_create_##type(void);			\
  int sorted_list_locate_##type(t_list_##type *list, type value);	\
  Int sorted_list_add_##type(t_list_##type *list, type value);		\


#define CREATE_SORTED_LIST_C(type)					\
  t_list_##type sorted_list_create_##type(void)				\
  {									\
    t_list_##type list;							\
    mmalloc(list.ptr, sizeof(type)*ALLOC_BLOCKSIZE);			\
    list.size = ALLOC_BLOCKSIZE;					\
    list.nelems = 0;							\
    return list;							\
  }									\
									\
  Int sorted_list_locate_##type(t_list_##type *list, type value)	\
  {									\
    Int l, u, m;							\
									\
    /* locate the range by bisection */					\
    l = 0;								\
    u = list->nelems;							\
									\
    while(u-l>64){							\
      m = (l+u)/2;							\
      if(list->ptr[m]>value){						\
	u=m;								\
      } else {								\
	l=m;								\
      }									\
    }									\
									\
    /* locate the value by linear search */				\
    while(l<u){								\
      if(list->ptr[l]>=value) break;					\
      l++;								\
    }									\
									\
    if(l<list->nelems && list->ptr[l]==value) return l;			\
    return list->nelems;						\
  }									\
  									\
  Int sorted_list_add_##type(t_list_##type *list, type value)		\
  {									\
									\
    Int l, retval;							\
    type temp;								\
    									\
    l = sorted_list_locate_##type(list, value);				\
    if(l<list->nelems && list->ptr[l]==value) return l;			\
									\
    /* check if we have enough of memory in the list */			\
    if(list->size==list->nelems-1){					\
      list->size += ALLOC_BLOCKSIZE;					\
      mrealloc(list->ptr, sizeof(type)*list->size,			\
	       ALLOC_BLOCKSIZE*sizeof(type));				\
    }									\
									\
    /* insert into array */						\
    retval = l;								\
    while(l<list->nelems){						\
      temp = list->ptr[l];						\
      list->ptr[l] = value;						\
      value = temp;							\
      l++;								\
    }									\
    list->ptr[list->nelems]=value;					\
    list->nelems++;							\
    return retval;							\
  }									\

