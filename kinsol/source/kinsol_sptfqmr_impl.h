/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:14:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPTFQMR linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSPTFQMR_IMPL_H
#define _KINSPTFQMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol_sptfqmr.h"

/*
 * -----------------------------------------------------------------
 * KINSptfqmr solver constant(s)
 * -----------------------------------------------------------------
 * KINSPTFQMR_MAXL : maximum dimension of Krylov subspace allowed by
 *                   default
 * -----------------------------------------------------------------
 */

#define KINSPTFQMR_MAXL 10

/*
 * -----------------------------------------------------------------
 * Types : struct KINSptfqmrMemRec and struct *KINSptfqmrMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *KINSptfqmrMem denotes a
 * pointer to a data structure of type struct KINSptfqmrMemRec. The
 * KINSptfqmrMemRec structure contains fields that must be accessed
 * by KINSPTFQMR/SPTFQMR solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct {

  /* problem specification data */

  int  q_maxl;          /* maximum allowable dimension of Krylov subspace      */
  int  q_pretype;       /* preconditioning type: PREC_NONE, PREC_RIGHT,
                           PREC_LEFT, or PREC_BOTH (used by SPTFQMR module and
                           defined in shared/include/iterative.h)              */
  booleantype q_new_uu; /* flag indicating if the iterate has been updated -
			   Jacobian must be updated/reevaluated (meant to be
			   used by user-supplied jtimes function)              */

  /* counters */

  long int q_nli;     /* number of linear iterations performed                */
  long int q_npe;     /* number of preconditioner evaluations                 */
  long int q_nps;     /* number of times preconditioner was applied to linear
		         system                                               */
  long int q_ncfl;    /* number of linear convergence failures                */
  long int q_nfeSG;   /* number of evaluations of the system function F(u) or
			 number of calls made to func routine                 */    
  long int q_njtimes; /* number of times the matrix-vector product J(u)*v
			 was computed or number of calls made to jtimes
			 routine                                              */

  /* functions (pointer) */

  KINSpilsPrecSetupFn q_pset;     /* routine called to compute preconditioner
				     matrix                                   */
  KINSpilsPrecSolveFn q_psolve;   /* subroutine called to solve a
				     preconditioned linear system             */ 
  KINSpilsJacTimesVecFn q_jtimes; /* function called to compute matrix-vector
				     product J(u)*v                           */

  /* memory references (pointers) */

  void *q_P_data; /* pointer to user-allocated memory block that is passed
		     to pset and psolve                                      */
  void *q_J_data; /* pointer to user-allocated memory block that is passed
		     to jtimes (only required if using a user-supplied
		     jtimes routine)                                         */

  SptfqmrMem q_sptfqmr_mem; /* pointer to SPTFQMR memory block (allocated by
			       SptfqmrMalloc routine)                        */

  /* miscellaneous data */

  int q_last_flag; /* last flag returned */

} KINSptfqmrMemRec, *KINSptfqmrMem;

/*
 * -----------------------------------------------------------------
 * KINSPTFQMR error messages
 * -----------------------------------------------------------------
 */

/* KINSptfqmr error messages */

#define KINSPTFQMR       "KINSptfqmr-- "
#define MSGQ_KINMEM_NULL KINSPTFQMR "KINSOL memory is NULL.\n\n"
#define MSGQ_MEM_FAIL    KINSPTFQMR "A memory request failed.\n\n"
#define MSGQ_BAD_NVECTOR KINSPTFQMR "A required vector operation is not implemented.\n\n"

/* KINSptfqmrSet* and KINSptfqmrGet* error messages */

#define KINSPTFQMR_SETGET       "KINSptfqmrSet*/KINSptfqmrGet*-- "
#define MSGQ_SETGET_KINMEM_NULL KINSPTFQMR_SETGET "KINSOL memory is NULL.\n\n"
#define MSGQ_SETGET_LMEM_NULL   KINSPTFQMR_SETGET "KINSPTFQMR memory is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
