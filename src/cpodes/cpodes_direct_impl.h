/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Common implementation header file for the CPDLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _CPDLS_IMPL_H
#define _CPDLS_IMPL_H

#include <cpodes/cpodes_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * CPDLS solver constants
 * -----------------------------------------------------------------
 * CPD_MSBJ   maximum number of steps between Jacobian evaluations
 * CPD_DGMAX  maximum change in gamma between Jacobian evaluations
 * -----------------------------------------------------------------
 */

#define CPD_MSBJ  50
#define CPD_DGMAX RCONST(0.2)

/*
 * -----------------------------------------------------------------
 * Types : CPDlsMemRec, CPDlsMem                             
 * -----------------------------------------------------------------
 * CPDlsMem is pointer to a CPDlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct {

  int d_type;             /* SUNDIALS_DENSE or SUNDIALS_BAND               */

  int d_n;                /* problem dimension                             */

  int d_ml;               /* lower bandwidth of Jacobian                   */
  int d_mu;               /* upper bandwidth of Jacobian                   */ 
  int d_smu;              /* upper bandwith of M = MIN(N-1,d_mu+d_ml)      */

  booleantype d_jacDQ;    /* SUNTRUE if using internal DQ Jacobian approx. */

  CPDlsDenseJacExplFn d_djacE; /* dense Jacobian routine (CP_EXPL)         */
  CPDlsDenseJacImplFn d_djacI; /* dense Jacobian routine (CP_IMPL)         */

  CPDlsBandJacExplFn d_bjacE;  /* band Jacobian routine (CP_EXPL)          */
  CPDlsBandJacImplFn d_bjacI;  /* band Jacobian routine (CP_IMPL)          */

  void *d_J_data;         /* data pointer passed to djac* or bjac*         */

  DlsMat d_M;             /* M = I - gamma * df/dy                         */
  DlsMat d_savedJ;        /* savedJ = old Jacobian                         */

  long *d_pivots;         /* pivots = pivot array for PM = LU              */
  
  long int d_nstlj;       /* nstlj = nst at last Jacobian eval.            */

  long int d_nje;         /* nje = no. of calls to jac                     */

  long int d_nfeDQ;       /* no. of calls to f due to DQ Jacobian approx.  */

  long int d_last_flag;   /* last error return flag                        */
  
} CPDlsMemRec, *CPDlsMem;

/*
 * -----------------------------------------------------------------
 * Types : CPDlsProjMemRec, CPDlsProjMem                             
 * -----------------------------------------------------------------
 * The type CPDlsProjMem is pointer to a CPDlsProjMemRec.
 * This structure contains CPDlsProj solver-specific data. 
 * -----------------------------------------------------------------
 */
  
typedef struct {

  int d_type;               /* always SUNDIALS_DENSE                         */

  int d_nc;                 /* number of constraints                         */
  int d_ny;                 /* number of states                              */

  booleantype d_jacPDQ;     /* SUNTRUE if using internal DQ Jacobian approx. */

  CPDlsDenseProjJacFn d_jacP; /* Jacobian routine to be called               */

  void *d_JP_data;          /* data pointer passed to jacP                   */

  int d_ftype;              /* factorization type (LU, QR, or SC)            */
  int d_pnorm;              /* projection norm (L2 or WRMS)                  */

  DlsMat d_G;               /* G = (dc/dy)^T, transpose of cnstr. Jacobian   */
  DlsMat d_savedG;          /* saved Jacobian (before factorization)         */

  int d_nr;                 /* no. of independent constraints (QRP)          */

  DlsMat d_K;               /* K matrix (s.p.d., form depends on ftype)      */
  long int *d_pivotsP;      /* pivotsP = pivot array (for ftype LU)          */

  realtype *d_beta;         /* beta array (for ftype QR)                     */

  realtype *d_wrk;          /* work array (for ftype QR or QRP)              */
  int d_len_wrk;            /* length of work array                          */

  long int d_nstljP;        /* nstljP = nst at last Jacobian eval.           */

  long int d_njeP;          /* njeP = no. of calls to jacP                   */

  long int d_nceDQ;         /* no. of calls to c due to DQ Jacobian approx.  */

} CPDlsProjMemRec, *CPDlsProjMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

int cpDlsDenseDQJacExpl(int N, realtype t,
			N_Vector y, N_Vector fy, 
			DlsMat Jac, void *data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
int cpDlsDenseDQJacImpl(int N, realtype t, realtype gm,
			N_Vector y, N_Vector yp, N_Vector r, 
			DlsMat Jac, void *data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int cpDlsBandDQJacExpl(int N, int mupper, int mlower,
		       realtype t, N_Vector y, N_Vector fy, 
		       DlsMat Jac, void *data,
		       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int cpDlsBandDQJacImpl(int N, int mupper, int mlower,
		       realtype t, realtype gm, 
		       N_Vector y, N_Vector yp, N_Vector r,
		       DlsMat Jac, void *data,
		       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int cpDlsDenseProjDQJac(int Nc, int Ny, realtype t,
			N_Vector y, N_Vector cy, 
			DlsMat Jac, void *data,
			N_Vector c_tmp1, N_Vector c_tmp2);

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGD_CPMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#define MSGD_BAD_FACT "fact_type has an illegal value."


#ifdef __cplusplus
}
#endif

#endif
