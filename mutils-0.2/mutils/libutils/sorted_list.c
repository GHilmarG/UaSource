/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "config.h"
#include "sorted_list.h"

void sorted_list_create(Int **list, Int *size)
{
  mmalloc((*list), sizeof(Int)*ALLOC_BLOCKSIZE);
  *size = ALLOC_BLOCKSIZE;
}

void sorted_list_create_double(Double **list, Int *size)
{
  mmalloc((*list), sizeof(Double)*ALLOC_BLOCKSIZE);
  *size = ALLOC_BLOCKSIZE;
}

void sorted_list_create_pair(Int **list, Double **listd, Int *size)
{
  mmalloc((*list) , sizeof(Int  )*ALLOC_BLOCKSIZE);
  mmalloc((*listd), sizeof(Double)*ALLOC_BLOCKSIZE);
  *size = ALLOC_BLOCKSIZE;
}

void sorted_list_add(Int **list, Int *nelems, Int *lsize, Int value)
{
  Int l, u;

  l = sorted_list_locate_Int(*list, *nelems, value);
  if(l<*nelems && (*list)[l]==value) return;

  /* check if we have enough of memory in the list */
  if(*nelems==*lsize){
    (*lsize) += ALLOC_BLOCKSIZE;
    mrealloc((*list), sizeof(Int)*(*lsize), sizeof(Int)*ALLOC_BLOCKSIZE);
  }
  
  /* insert into array */
  u = *nelems;
  while(l!=u){
    (*list)[u] = (*list)[u-1];
    u--;
  }
  
  (*list)[l] = value;
  (*nelems)++;
}

void sorted_list_add_accum(Int **list, Int *nelems, Int *lsize, Int value, Double **dlist, Double dvalue)
{
  Int l, temp;
  Double dtemp;

  l = sorted_list_locate_Int(*list, *nelems, value);
  if(l<*nelems && (*list)[l]==value) {
    (*dlist)[l]+=dvalue;
    return;
  }

  /* check if we have enough of memory in the list */
  if(*nelems==*lsize){
    (*lsize) += ALLOC_BLOCKSIZE;
    mrealloc((*list), sizeof(Int)*(*lsize), ALLOC_BLOCKSIZE*sizeof(Int));
    mrealloc((*dlist), sizeof(Double)*(*lsize), ALLOC_BLOCKSIZE*sizeof(Double));
  }

  /* insert into array */
  while(l<*nelems){

    temp = (*list)[l];
    (*list)[l] = value;
    value = temp;

    dtemp = (*dlist)[l];
    (*dlist)[l] = dvalue;
    dvalue = dtemp;

    l++;
  }
  (*list)[(*nelems)]=value;
  (*dlist)[(*nelems)]=dvalue;
  (*nelems)++;
}


/* define the inline functions here if we do not have INLINE */
#if !HAVE_INLINE

#define C DEFINITION OF SORTED LIST

SORTED_LIST_LOCATE_C(Int);
SORTED_LIST_LOCATE_C(Long);
SORTED_LIST_LOCATE_C(dimType);
SORTED_LIST_LOCATE_C(indexType);

SORTED_LIST_ADD_STATIC_C(Int);
SORTED_LIST_ADD_STATIC_C(dimType);

SORTED_LIST_ADD_STATIC_ACCUM_C(Int, Double);
SORTED_LIST_ADD_STATIC_ACCUM_C(dimType, Double);

#endif


