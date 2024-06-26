/* ========================================================================== */
/* === UMF_valid_symbolic =================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"
#include "umf_valid_symbolic.h"

/* Returns TRUE if the Symbolic object is valid, FALSE otherwise. */
/* The UMFPACK_report_symbolic routine does a more thorough check. */

GLOBAL Int UMF_valid_symbolic
(
    SymbolicType *Symbolic
)
{
    /* This routine does not check the contents of the individual arrays, so */
    /* it can miss some errors.  All it checks for is the presence of the */
    /* arrays, and the Symbolic "valid" entry. */

    if (!Symbolic)
    {
	return (FALSE) ;
    }

    if (Symbolic->valid != SYMBOLIC_VALID)
    {
	/* Symbolic does not point to a SymbolicType object */
	return (FALSE) ;
    }

    if (!Symbolic->Cperm_init || !Symbolic->Rperm_init ||
	!Symbolic->Front_npivcol || !Symbolic->Front_1strow ||
	!Symbolic->Front_leftmostdesc ||
	!Symbolic->Front_parent || !Symbolic->Chain_start ||
	!Symbolic->Chain_maxrows || !Symbolic->Chain_maxcols ||
	Symbolic->n_row <= 0 || Symbolic->n_col <= 0)
    {
	return (FALSE) ;
    }

    return (TRUE) ;
}
