/* ========================================================================== */
/* === UMFPACK_free_numeric ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

/*  User-callable.  Free the entire Numeric object (consists of 11 to 13
 *  malloc'd objects.  See UMFPACK_free_numeric.h for details.
 */

#include "umf_internal.h"
#include "umf_free.h"

GLOBAL void UMFPACK_free_numeric
(
    void **NumericHandle
)
{

    NumericType *Numeric ;
    if (!NumericHandle)
    {
	return ;
    }
    Numeric = *((NumericType **) NumericHandle) ;
    if (!Numeric)
    {
	return ;
    }

    /* these 9 objects always exist */
    (void) UMF_free ((void *) Numeric->D) ;
    (void) UMF_free ((void *) Numeric->Rperm) ;
    (void) UMF_free ((void *) Numeric->Cperm) ;
    (void) UMF_free ((void *) Numeric->Lpos) ;
    (void) UMF_free ((void *) Numeric->Lilen) ;
    (void) UMF_free ((void *) Numeric->Lip) ;
    (void) UMF_free ((void *) Numeric->Upos) ;
    (void) UMF_free ((void *) Numeric->Uilen) ;
    (void) UMF_free ((void *) Numeric->Uip) ;

    /* Rs does not exist if scaling was not performed */
    (void) UMF_free ((void *) Numeric->Rs) ;

    /* Upattern can only exist for singular or rectangular matrices */
    (void) UMF_free ((void *) Numeric->Upattern) ;

    /* these 2 objects always exist */
    (void) UMF_free ((void *) Numeric->Memory) ;
    (void) UMF_free ((void *) Numeric) ;

    *NumericHandle = (void *) NULL ;
}
