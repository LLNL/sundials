/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * Common implementation header file for the ARKDLS linear solvers.
 *--------------------------------------------------------------*/

#ifndef _ARKDLS_IMPL_H
#define _ARKDLS_IMPL_H

#include <arkode/arkode_direct.h>

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

  int d_type;             /* SUNDIALS_DENSE or SUNDIALS_BAND              */

  long int d_n;           /* problem dimension                            */

  long int d_ml;          /* lower bandwidth of Jacobian                  */
  long int d_mu;          /* upper bandwidth of Jacobian                  */ 
  long int d_smu;         /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */

  booleantype d_jacDQ;    /* TRUE if using internal DQ Jacobian approx.   */
  ARKDlsDenseJacFn d_djac; /* dense Jacobian routine to be called          */
  ARKDlsBandJacFn d_bjac;  /* band Jacobian routine to be called           */
  void *d_J_data;         /* user data is passed to djac or bjac          */

  DlsMat d_M;             /* M = I - gamma * df/dy                        */
  DlsMat d_savedJ;        /* savedJ = old Jacobian                        */

  int *d_pivots;          /* pivots = int pivot array for PM = LU         */
  long int *d_lpivots;    /* lpivots = long int pivot array for PM = LU   */

  long int  d_nstlj;      /* nstlj = nst at last Jacobian eval.           */

  long int d_nje;         /* nje = no. of calls to jac                    */

  long int d_nfeDQ;       /* no. of calls to f due to DQ Jacobian approx. */

  long int d_last_flag;   /* last error return flag                       */
  
} *ARKDlsMem;


/*---------------------------------------------------------------
 Types: ARKDlsMassMemRec, ARKDlsMassMem

 ARKDlsMassMem is pointer to a ARKDlsMassMemRec structure.
---------------------------------------------------------------*/
typedef struct ARKDlsMassMemRec {

  int d_type;                /* SUNDIALS_DENSE or SUNDIALS_BAND            */

  long int d_n;              /* problem dimension                          */

  long int d_ml;             /* lower bandwidth of mass matrix             */
  long int d_mu;             /* upper bandwidth of mass matrix             */ 
  long int d_smu;            /* upper bandwith of M = MIN(N-1,d_mu+d_ml)   */

  ARKDlsDenseMassFn d_dmass; /* dense mass matrix routine to be called     */
  ARKDlsBandMassFn d_bmass;  /* band mass matrix routine to be called      */
  void *d_M_data;            /* user data is passed to djac or bjac        */

  DlsMat d_M;                /* mass matrix structure                      */
  DlsMat d_M_lu;             /* mass matrix structure for LU decomposition */

  int *d_pivots;             /* pivots = int pivot array for PM = LU       */
  long int *d_lpivots;       /* lpivots = long int pivot array for PM = LU */

  long int d_nme;            /* nje = no. of calls to mass matrix routine  */

  long int d_last_flag;      /* last error return flag                     */
  
} *ARKDlsMassMem;


/*---------------------------------------------------------------
 Prototypes of internal functions
---------------------------------------------------------------*/
int arkDlsDenseDQJac(long int N, realtype t, N_Vector y, 
                     N_Vector fy, DlsMat Jac, void *data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int arkDlsBandDQJac(long int N, long int mupper, long int mlower,
                    realtype t, N_Vector y, N_Vector fy, 
                    DlsMat Jac, void *data, N_Vector tmp1, 
                    N_Vector tmp2, N_Vector tmp3);

/* Auxilliary functions */

int arkDlsInitializeCounters(ARKDlsMem arkdls_mem);


/*---------------------------------------------------------------
 Error Messages
---------------------------------------------------------------*/
#define MSGD_ARKMEM_NULL    "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSGD_BAD_SIZES      "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL       "A memory request failed."
#define MSGD_LMEM_NULL      "Linear solver memory is NULL."
#define MSGD_MASSMEM_NULL   "Mass matrix solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGD_MASSFUNC_FAILED "The mass matrix routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
