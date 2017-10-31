/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Radu Serban @ LLNL
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
 * Common implementation header file for the KINDLS linear solvers.
 * -----------------------------------------------------------------*/

#ifndef _KINDLS_IMPL_H
#define _KINDLS_IMPL_H

#include <kinsol/kinsol_direct.h>
#include "kinsol_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*------------------------------------------------------------------
  Types: KINDlsMemRec, KINDlsMem                             
  ------------------------------------------------------------------
  The type KINDlsMem is pointer to a KINDlsMemRec.
  This structure contains KINDLS solver-specific data. 
  ------------------------------------------------------------------*/

typedef struct KINDlsMemRec {

  booleantype jacDQ;   /* SUNTRUE if using internal DQ Jacobian approx. */
  KINDlsJacFn jac;     /* Jacobian routine to be called                 */
  void *J_data;        /* J_data is passed to jac                       */
  
  SUNLinearSolver LS;  /* generic direct linear solver object           */

  SUNMatrix J;         /* problem Jacobian                              */

  long int nje  ;      /* no. of calls to jac                           */
    
  long int nfeDQ;      /* no. of calls to F due to DQ Jacobian approx.  */
    
  long int last_flag;  /* last error return flag                        */
    
} *KINDlsMem;


/*------------------------------------------------------------------
  Prototypes of internal functions
  ------------------------------------------------------------------*/

/* difference-quotient Jacobian approximation routines */
int kinDlsDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac,
                void *data, N_Vector tmp1, N_Vector tmp2);

int kinDlsDenseDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac,
                     KINMem kin_mem, N_Vector tmp1, N_Vector tmp2);

int kinDlsBandDQJac(N_Vector u, N_Vector fu, SUNMatrix Jac,
                    KINMem kin_mem, N_Vector tmp1, N_Vector tmp2);

/* generic linit/lsetup/lsolve/lfree interface routines for KINSOL to call */
int kinDlsInitialize(KINMem kin_mem);

int kinDlsSetup(KINMem kin_mem);

int kinDlsSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                realtype *sJpnorm, realtype *sFdotJp);
                  
int kinDlsFree(KINMem kin_mem);

/* Auxilliary functions */
int kinDlsInitializeCounters(KINDlsMem kindls_mem);

/*------------------------------------------------------------------
  Error Messages
  ------------------------------------------------------------------*/

#define MSGD_KINMEM_NULL "KINSOL memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_MEM_FAIL    "A memory request failed."
#define MSGD_LMEM_NULL   "Linear solver memory is NULL."
#define MSGD_BAD_SIZES   "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGD_MATZERO_FAILED "The SUNMatZero routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
