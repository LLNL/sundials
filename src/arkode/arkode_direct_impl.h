/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Implementation header file for the ARKDLS linear solver 
 * interface
 *--------------------------------------------------------------*/

#ifndef _ARKDLS_IMPL_H
#define _ARKDLS_IMPL_H

#include <arkode/arkode_direct.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
 ARKDLS solver constants:

 ARKD_MSBJ   maximum number of steps between Jacobian evaluations
 ARKD_DGMAX  maximum change in gamma between Jacobian evaluations
---------------------------------------------------------------*/
#define ARKD_MSBJ  50
#define ARKD_DGMAX RCONST(0.2)


/*---------------------------------------------------------------
 Types: ARKDlsMemRec, ARKDlsMem                             

 ARKDlsMem is pointer to a ARKDlsMemRec structure.
---------------------------------------------------------------*/
typedef struct ARKDlsMemRec {

  booleantype jacDQ;    /* SUNTRUE if using internal DQ Jacobian approx. */
  ARKDlsJacFn jac;      /* Jacobian routine to be called                 */
  void *J_data;         /* user data is passed to jac                    */

  SUNLinearSolver LS;   /* generic direct linear solver object           */

  SUNMatrix A;          /* A = M - gamma * df/dy                         */
  SUNMatrix savedJ;     /* savedJ = old Jacobian                         */

  N_Vector x;           /* solution vector used by SUNLinearSolver       */
  
  long int nstlj;       /* nstlj = nst at last Jacobian eval.            */

  long int nje;         /* nje = no. of calls to jac                     */

  long int nfeDQ;       /* no. of calls to f due to DQ Jacobian approx.  */

  long int last_flag;   /* last error return flag                        */
  
} *ARKDlsMem;


/*---------------------------------------------------------------
 Types: ARKDlsMassMemRec, ARKDlsMassMem

 ARKDlsMassMem is pointer to a ARKDlsMassMemRec structure.
---------------------------------------------------------------*/
typedef struct ARKDlsMassMemRec {

  ARKDlsMassFn mass;      /* user-provided mass matrix routine to call  */

  SUNLinearSolver LS;     /* generic direct linear solver object        */

  SUNMatrix M;            /* mass matrix structure                      */
  SUNMatrix M_lu;         /* mass matrix structure for LU decomposition */

  N_Vector x;             /* solution vector used by SUNLinearSolver    */

  booleantype time_dependent;  /* flag stating whether M depends on t   */
  
  long int mass_setups;   /* number of mass matrix-solver setup calls   */

  long int mass_solves;   /* number of mass matrix solve calls          */

  long int mass_mults;    /* number of mass matrix product calls        */

  long int last_flag;     /* last error return flag                     */
  
} *ARKDlsMassMem;


/*---------------------------------------------------------------
 Prototypes of internal functions
---------------------------------------------------------------*/

/* difference-quotient Jacobian approximation routines */
int arkDlsDQJac(realtype t, N_Vector y, N_Vector fy, 
                SUNMatrix Jac, void *data, N_Vector tmp1, 
                N_Vector tmp2, N_Vector tmp3);
int arkDlsDenseDQJac(realtype t, N_Vector y, N_Vector fy, 
                     SUNMatrix Jac, ARKodeMem ark_mem, N_Vector tmp1);
int arkDlsBandDQJac(realtype t, N_Vector y, N_Vector fy, 
                    SUNMatrix Jac, ARKodeMem ark_mem, N_Vector tmp1, 
                    N_Vector tmp2);

/* generic linit/lsetup/lsolve/lfree interface routines for ARKode to call */
int arkDlsInitialize(ARKodeMem ark_mem);

int arkDlsSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                N_Vector fpred, booleantype *jcurPtr, 
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

int arkDlsSolve(ARKodeMem ark_mem, N_Vector b, N_Vector ycur, N_Vector fcur);

int arkDlsFree(ARKodeMem ark_mem);

/* generic minit/msetup/mmult/msolve/mfree routines for ARKode to call */  
int arkDlsMassInitialize(ARKodeMem ark_mem);
  
int arkDlsMassSetup(ARKodeMem ark_mem, N_Vector vtemp1,
                    N_Vector vtemp2, N_Vector vtemp3); 

int arkDlsMassMult(ARKodeMem ark_mem, N_Vector v, N_Vector Mv);

int arkDlsMassSolve(ARKodeMem ark_mem, N_Vector b);

int arkDlsMassFree(ARKodeMem ark_mem);

/* Auxilliary functions */
int arkDlsInitializeCounters(ARKDlsMem arkdls_mem);

int arkDlsInitializeMassCounters(ARKDlsMassMem arkdls_mem);


/*---------------------------------------------------------------
 Error Messages
---------------------------------------------------------------*/
#define MSGD_ARKMEM_NULL         "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR         "A required vector operation is not implemented."
#define MSGD_BAD_SIZES           "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL            "A memory request failed."
#define MSGD_LMEM_NULL           "Linear solver memory is NULL."
#define MSGD_MASSMEM_NULL        "Mass matrix solver memory is NULL."
#define MSGD_JACFUNC_FAILED      "The Jacobian routine failed in an unrecoverable manner."
#define MSGD_MASSFUNC_FAILED     "The mass matrix routine failed in an unrecoverable manner."
#define MSGD_MATCOPY_FAILED      "The SUNMatCopy routine failed in an unrecoverable manner."
#define MSGD_MATZERO_FAILED      "The SUNMatZero routine failed in an unrecoverable manner."
#define MSGD_MATSCALEADD_FAILED  "The SUNMatScaleAdd routine failed in an unrecoverable manner."
#define MSGD_MATSCALEADDI_FAILED "The SUNMatScaleAddI routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
