/*
 * -----------------------------------------------------------------
 * $Revision: 1.5.2.1 $
 * $Date: 2005-01-26 22:05:17 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPGMR linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSPGMR_IMPL_H
#define _KINSPGMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinspgmr.h"
#include "nvector.h"
#include "iterative.h"
#include "spgmr.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * KINSpgmr solver constant
 * -----------------------------------------------------------------
 * KINSPGMR_MAXL : maximum dimension of Krylov subspace allowed by
 *                 default
 * -----------------------------------------------------------------
 */

#define KINSPGMR_MAXL 10

/*
 * -----------------------------------------------------------------
 * Types : struct KINSpgmrMemRec and struct *KINSpgmrMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *KINSpgmrMem denotes a
 * pointer to a data structure of type struct KINSpgmrMemRec. The
 * KINSpgmrMemRec structure contains fields that must be accessible
 * by KINSPGMR/SPGMR solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct {

  /* problem specification data */

  int  g_maxl;          /* maximum allowable dimension of Krylov subspace      */     
  int  g_pretype;       /* preconditioning type: PREC_NONE, PREC_RIGHT,
			   PREC_LEFT or PREC_BOTH (used by SPGMR module and
			   defined in shared/include/iterative.h)              */
  int  g_gstype;        /* Gram-Schmidt orthogonalization procedure:
			   CLASSICAL_GS or MODIFIED_GS (used by SPGMR module
			   and defined in shared/include/iterative.h)          */
  booleantype g_new_uu; /* flag indicating if the iterate has been updated -
			   Jacobian must be updated/reevaluated (meant to be
			   used by user-supplied jtimes function)              */
  int g_maxlrst;        /* maximum number of times the SPGMR linear solver
			   can be restarted                                    */

  /* counters */

  long int g_nli;     /* number of linear iterations performed                 */
  long int g_npe;     /* number of preconditioner evaluations                  */
  long int g_nps;     /* number of times preconditioner was applied to linear
		         system                                                */
  long int g_ncfl;    /* number of linear convergence failures                 */
  long int g_nfeSG;   /* number of evaluations of the system function F(u) or
			 number of calls made to func routine                  */    
  long int g_njtimes; /* number of times the matrix-vector product J(u)*v
			 was computed or number of calls made to jtimes
			 routine                                               */

  /* functions */

  KINSpgmrPrecSetupFn g_pset;     /* routine called to compute preconditioner
				     matrix                                    */
  KINSpgmrPrecSolveFn g_psolve;   /* subroutine called to solve a
				     preconditioned linear system              */ 
  KINSpgmrJacTimesVecFn g_jtimes; /* function called to compute matrix-vector
				     product J(u)*v                            */

  /* memory references (pointers) */

  void *g_P_data; /* pointer to user-allocated memory block that is passed
		     to pset and psolve                                        */
  void *g_J_data; /* pointer to user-allocated memory block that is passed
		     to jtimes (only required if using a user-supplied
		     jtimes routine)                                           */

  SpgmrMem g_spgmr_mem; /* pointer to SPGMR memory block (allocated by
			   SpgmrMalloc routine)                                */

  /* miscellaneous data */

  int g_last_flag; /* last flag returned                                       */

} KINSpgmrMemRec, *KINSpgmrMem;

/*
 * -----------------------------------------------------------------
 * KINSPGMR error messages
 * -----------------------------------------------------------------
 */

/* KINSpgmr error messages */

#define KINSPGMR        "KINSpgmr-- "
#define MSGS_KINMEM_NULL KINSPGMR "KINSOL memory is NULL.\n\n"
#define MSGS_MEM_FAIL    KINSPGMR "A memory request failed.\n\n"
#define MSGS_BAD_NVECTOR KINSPGMR "A required vector operation is not implemented.\n\n"

/* KINSpgmrSet* and KINSpgmrGet* error messages */

#define KINSPGMR_SETGET        "KINSpgmrSet*/KINSpgmrGet*-- "
#define MSGS_SETGET_KINMEM_NULL KINSPGMR_SETGET "KINSOL memory is NULL. \n\n"
#define MSGS_SETGET_LMEM_NULL   KINSPGMR_SETGET "KINSPGMR memory is NULL.\n\n"

/* KINSpgmrSetMaxRestarts error message */

#define MSGS_KINS_NEG_MAXRS "KINSpgmrSetMaxRestarts-- maxrs < 0 illegal.\n\n"

#ifdef __cplusplus
}
#endif

#endif
