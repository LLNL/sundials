/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Common implementation header file for the scaled, preconditioned
 * linear solver modules.
 * -----------------------------------------------------------------*/

#ifndef _KINSPILS_IMPL_H
#define _KINSPILS_IMPL_H

#include <kinsol/kinsol_spils.h>
#include "kinsol_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*------------------------------------------------------------------
  keys for KINPrintInfo (do not use 1 -> conflict with PRNT_RETVAL)
  ------------------------------------------------------------------*/

#define PRNT_NLI   101
#define PRNT_EPS   102


/*------------------------------------------------------------------
  Types : struct KINSpilsMemRec, struct *KINSpilsMem
  ------------------------------------------------------------------
  A variable declaration of type struct *KINSpilsMem denotes a
  pointer to a data structure of type struct KINSpilsMemRec. The
  KINSpilsMemRec structure contains fields that must be accessible
  by SPILS module routines.
  ------------------------------------------------------------------*/

typedef struct KINSpilsMemRec {

  long int npe;       /* npe = total number of precond calls          */
  long int nli;       /* nli = total number of linear iterations      */
  long int nps;       /* nps = total number of psolve calls           */
  long int ncfl;      /* ncfl = total number of convergence failures  */
  long int nfes;      /* nres = total number of calls to F(u)         */
  long int njtimes;   /* njtimes = total number of calls to jtimes    */

  booleantype new_uu; /* flag indicating if the iterate has been 
                         updated - the Jacobian must be updated or 
                         reevaluated (meant to be used by a
                         user-supplied jtimes function                */

  SUNLinearSolver LS; /* generic iterative linear solver object       */

  long int last_flag; /* last error return flag                       */

  /* Preconditioner computation
     (a) user-provided:
         - pdata == user_data
         - pfree == NULL (the user dealocates memory)
     (b) internal preconditioner module
         - pdata == kin_mem
         - pfree == set by the prec. module and called in kinSpilsFree */
  KINSpilsPrecSetupFn pset;
  KINSpilsPrecSolveFn psolve;
  int (*pfree)(KINMem kin_mem);
  void *pdata;

  /* Jacobian times vector compuation
     (a) jtimes function provided by the user:
         - jdata == user_data
         - jtimesDQ == SUNFALSE
     (b) internal jtimes
         - jdata == kin_mem
         - jtimesDQ == SUNTRUE */
  booleantype jtimesDQ;
  KINSpilsJacTimesVecFn jtimes;
  void *jdata;

} *KINSpilsMem;


/*------------------------------------------------------------------
  Prototypes of internal functions
  ------------------------------------------------------------------*/

/* Interface routines called by system SUNLinearSolvers */
int KINSpilsATimes(void *kin_mem, N_Vector v, N_Vector z);
int KINSpilsPSetup(void *kin_mem);
int KINSpilsPSolve(void *kin_mem, N_Vector r, N_Vector z,
                   realtype tol, int lr);

/* Difference quotient approximation for Jacobian times vector */
int KINSpilsDQJtimes(N_Vector v, N_Vector Jv,
                     N_Vector u, booleantype *new_u,
                     void *data);

/* Generic linit/lsetup/lsolve/lfree interface routines for KINSOL to call */
int kinSpilsInitialize(KINMem kin_mem);

int kinSpilsSetup(KINMem kin_mem);

int kinSpilsSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                  realtype *sJpnorm, realtype *sFdotJp);

int kinSpilsFree(KINMem kin_mem);

/* Auxilliary functions */
int kinSpilsInitializeCounters(KINSpilsMem kinspils_mem);


/*------------------------------------------------------------------
  Error messages
  ------------------------------------------------------------------*/

#define MSGS_KINMEM_NULL "KINSOL memory is NULL."
#define MSGS_MEM_FAIL    "A memory request failed."
#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGS_LMEM_NULL   "Linear solver memory is NULL."
#define MSGS_NEG_MAXRS   "maxrs < 0 illegal."

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."


/*------------------------------------------------------------------
  Info messages
  ------------------------------------------------------------------*/

#define INFO_NLI  "nli_inc = %d"

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define INFO_EPS  "residual norm = %12.3Lg  eps = %12.3Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define INFO_EPS  "residual norm = %12.3lg  eps = %12.3lg"

#else

#define INFO_EPS  "residual norm = %12.3g  eps = %12.3g"

#endif


#ifdef __cplusplus
}
#endif

#endif
