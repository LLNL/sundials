/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-01-24 23:55:50 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPBCG linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSPBCG_IMPL_H
#define _KINSPBCG_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol.h"
#include "kinspbcg.h"
#include "nvector.h"
#include "spbcg.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * KINSpbcg solver constant(s)
 * -----------------------------------------------------------------
 * KINSPBCG_MAXL : maximum dimension of Krylov subspace allowed by
 *                 default
 * -----------------------------------------------------------------
 */

#define KINSPBCG_MAXL 10

/*
 * -----------------------------------------------------------------
 * Types : struct KINSpbcgMemRec and struct *KINSpbcgMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *KINSpbcgMem denotes a
 * pointer to a data structure of type struct KINSpbcgMemRec. The
 * KINSpbcgMemRec structure contains fields that must be accessed
 * by KINSPBCG/SPBCG solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct {

  /* problem specification data */

  int  b_maxl;          /* maximum allowable dimension of Krylov subspace    */     
  int  b_pretype;       /* preconditioning type: PREC_NONE, PREC_RIGHT,
                           PREC_LEFT, or PREC_BOTH (used by SPBCG module and
                           defined in shared/include/iterative.h)            */
  booleantype b_new_uu; /* flag indicating if the iterate has been updated -
			   Jacobian must be updated/reevaluated (meant to be
			   used by user-supplied jtimes function)            */

  /* counters */

  long int b_nli;     /* number of linear iterations performed                */
  long int b_npe;     /* number of preconditioner evaluations                 */
  long int b_nps;     /* number of times preconditioner was applied to linear
		         system                                               */
  long int b_ncfl;    /* number of linear convergence failures                */
  long int b_nfeSG;   /* number of evaluations of the system function F(u) or
			 number of calls made to func routine                 */    
  long int b_njtimes; /* number of times the matrix-vector product J(u)*v
			 was computed or number of calls made to jtimes
			 routine                                              */

  /* functions (pointer) */

  KINSpbcgPrecSetupFn b_pset;     /* routine called to compute preconditioner
				     matrix                                   */
  KINSpbcgPrecSolveFn b_psolve;   /* subroutine called to solve a
				     preconditioned linear system             */ 
  KINSpbcgJacTimesVecFn b_jtimes; /* function called to compute matrix-vector
				     product J(u)*v                           */

  /* memory references (pointers) */

  void *b_P_data; /* pointer to user-allocated memory block that is passed
		     to pset and psolve                                    */
  void *b_J_data; /* pointer to user-allocated memory block that is passed
		     to jtimes (only required if using a user-supplied
		     jtimes routine)                                       */

  SpbcgMem b_spbcg_mem; /* pointer to SPBCG memory block (allocated by
			   SpbcgMalloc routine)                            */

  /* miscellaneous data */

  int b_last_flag; /* last flag returned */

} KINSpbcgMemRec, *KINSpbcgMem;

/*
 * -----------------------------------------------------------------
 * KINSPBCG error messages
 * -----------------------------------------------------------------
 */

/* KINSpbcg error messages */

#define KINSPBCG         "KINSpbcg-- "
#define MSGB_KINMEM_NULL KINSPBCG "KINSOL memory is NULL.\n\n"
#define MSGB_MEM_FAIL    KINSPBCG "A memory request failed.\n\n"
#define MSGB_BAD_NVECTOR KINSPBCG "A required vector operation is not implemented.\n\n"

/* KINSpbcgSet* and KINSpbcgGet* error messages */

#define KINSPBCG_SETGET         "KINSpbcgSet*/KINSpbcgGet*-- "
#define MSGB_SETGET_KINMEM_NULL KINSPBCG_SETGET "KINSOL memory is NULL.\n\n"
#define MSGB_SETGET_LMEM_NULL   KINSPBCG_SETGET "KINSPBCG memory is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
