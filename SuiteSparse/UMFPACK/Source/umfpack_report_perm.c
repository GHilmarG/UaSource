/* ========================================================================== */
/* === UMFPACK_report_perm ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints a permutation vector.  See umfpack_report_perm.h
    for details.

    Dynamic memory usage:  Allocates a size max(np,1)*sizeof(Int) workspace via
    a single call to UMF_malloc and then frees all of it via UMF_free on return.
*/

#include "umf_internal.h"
#include "umf_report_perm.h"
#include "umf_malloc.h"
#include "umf_free.h"

GLOBAL Int UMFPACK_report_perm
(
    Int np,
    const Int Perm [ ],
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl, *W, status ;

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    W = (Int *) UMF_malloc (MAX (np,1), sizeof (Int)) ;
    status = UMF_report_perm (np, Perm, W, prl, 1) ;
    (void) UMF_free ((void *) W) ;
    return (status) ;
}
