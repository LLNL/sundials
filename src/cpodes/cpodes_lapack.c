/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:06 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a CPODES dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/*
 * NOTE: the only operations (all O(n)) that do not use Blas/Lapack 
 *       functions are:
 *
 *   - matrix plus identity (I-gamma*J)
 *     (in lsetup)
 *   - diagonal matrix times vector (D^(-1)*x) 
 *     (in lsolveP for QR, QRP, and SC)
 *   - permutation of a diagonal matrix (P*D*P^T) 
 *     (in lsolveP for LU; P uses pivots from dgetrv)
 *   - permutation matrix times vector (P^T*x)
 *     (in lsolveP for LU; P uses pivots from dgetr)
 *   - permutation matrix times vector (P^T*b)
 *     (in lsolveP for QRP; P uses pivots from dgeqp3)
 */

/*
 * TODO:
 *   
 *   cplLUcomputeKD
 *   cplQRcomputeKD
 *   cplSCcomputeKD
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_lapack_impl.h"
#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* CPLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int cpLapackDenseInit(CPodeMem cp_mem);
static int cpLapackDenseSetup(CPodeMem cp_mem, int convfail, 
                              N_Vector yP, N_Vector ypP, N_Vector fctP, 
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpLapackDenseFree(CPodeMem cp_mem);
static void cpLapackDenseAddI(LapackMat A);

/* CPLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int cpLapackBandInit(CPodeMem cp_mem);
static int cpLapackBandSetup(CPodeMem cp_mem, int convfail, 
                             N_Vector yP, N_Vector ypP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpLapackBandFree(CPodeMem cp_mem);
static void cpLapackBandAddI(LapackMat A);


/* CPLAPACK DENSE linitP, lsetupP, lsolveP, lmultP, and lfreeP routines */
static int cpLapackDenseProjInit(CPodeMem cp_mem);
static int cpLapackDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1);
static int cpLapackDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                                  N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector s_tmp1);
static void cpLapackDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx);
static void cpLapackDenseProjFree(CPodeMem cp_mem);

/* Private functions for LU, QR, and SC projection */
static void cplLUcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cplQRcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cplSCcomputeKD(CPodeMem cp_mem, N_Vector d);

/* CPLAPACK DENSE DQ integration Jacobian functions */
static int cpLapackDenseDQJacExpl(int N, realtype t,
                                  N_Vector y, N_Vector fy, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackDenseDQJacImpl(int N, realtype t, realtype gm,
                                  N_Vector y, N_Vector yp, N_Vector r, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* CPLAPACK BAND DQ integration Jacobian functions */
static int cpLapackBandDQJacExpl(int N, int mupper, int mlower,
                                 realtype t, N_Vector y, N_Vector fy, 
                                 LapackMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpLapackBandDQJacImpl(int N, int mupper, int mlower,
                                 realtype t, realtype gm, 
                                 N_Vector y, N_Vector yp, N_Vector r,
                                 LapackMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* CPLAPACK DQ projection Jacobian functions */
static int cpLapackDenseProjDQJac(int Nc, int Ny, realtype t,
                                  N_Vector y, N_Vector cy, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector c_tmp1, N_Vector c_tmp2);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type       (cp_mem->cp_ode_type)
#define lmm_type       (cp_mem->cp_lmm_type)
#define fe             (cp_mem->cp_fe)
#define fi             (cp_mem->cp_fi)
#define f_data         (cp_mem->cp_f_data)
#define uround         (cp_mem->cp_uround)
#define nst            (cp_mem->cp_nst)
#define tn             (cp_mem->cp_tn)
#define h              (cp_mem->cp_h)
#define gamma          (cp_mem->cp_gamma)
#define gammap         (cp_mem->cp_gammap)
#define gamrat         (cp_mem->cp_gamrat)
#define ewt            (cp_mem->cp_ewt)

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lfree          (cp_mem->cp_lfree)
#define lmem           (cp_mem->cp_lmem)
#define tempv          (cp_mem->cp_tempv)
#define lsetup_exists  (cp_mem->cp_lsetup_exists)

#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)
#define pnorm          (cp_mem->cp_proj_norm)

#define linitP         (cp_mem->cp_linitP)
#define lsetupP        (cp_mem->cp_lsetupP)
#define lsolveP        (cp_mem->cp_lsolveP)
#define lmultP         (cp_mem->cp_lmultP)
#define lfreeP         (cp_mem->cp_lfreeP)
#define lmemP          (cp_mem->cp_lmemP)
#define lsetupP_exists (cp_mem->cp_lsetupP_exists)

#define mtype          (cplapack_mem->l_mtype)
#define n              (cplapack_mem->l_n)
#define ml             (cplapack_mem->b_ml)
#define mu             (cplapack_mem->b_mu)
#define smu            (cplapack_mem->b_smu)
#define djacE          (cplapack_mem->d_jacE)
#define djacI          (cplapack_mem->d_jacI)
#define bjacE          (cplapack_mem->b_jacE)
#define bjacI          (cplapack_mem->b_jacI)
#define M              (cplapack_mem->l_M)
#define savedJ         (cplapack_mem->l_savedJ)
#define pivots         (cplapack_mem->l_pivots)
#define nstlj          (cplapack_mem->l_nstlj)
#define nje            (cplapack_mem->l_nje)
#define nfeDQ          (cplapack_mem->l_nfeDQ)
#define J_data         (cplapack_mem->l_J_data)
#define last_flag      (cplapack_mem->l_last_flag)

#define ny             (cplapackP_mem->l_ny)
#define nc             (cplapackP_mem->l_nc)
#define nr             (cplapackP_mem->l_nr)
#define djacP          (cplapackP_mem->l_jacP)
#define JP_data        (cplapackP_mem->l_JP_data)
#define ftype          (cplapackP_mem->l_ftype)
#define G              (cplapackP_mem->l_G)
#define savedG         (cplapackP_mem->l_savedG)
#define K              (cplapackP_mem->l_K)
#define pivotsP        (cplapackP_mem->l_pivotsP)
#define beta           (cplapackP_mem->l_beta)
#define wrk            (cplapackP_mem->l_wrk)
#define len_wrk        (cplapackP_mem->l_len_wrk)
#define nstljP         (cplapackP_mem->l_nstljP)
#define njeP           (cplapackP_mem->l_njeP)
#define nceDQ          (cplapackP_mem->l_nceDQ)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CPLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  CPLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cp_linit, cp_lsetup, cp_lsolve, cp_lfree fields in (*cpode_mem)
 * to be cpLapackDenseInit, cpLapackDenseSetup, cpLapackDenseSolve, 
 * and cpLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type CPLapackMemRec and sets the cp_lmem field in 
 * (*cpode_mem) to the address of this structure.  It sets lsetup_exists 
 * in (*cpode_mem) to TRUE, and the d_jac field to the default 
 * cpLapackDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int CPLapackDense(void *cpode_mem, int N)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackDense", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPLAPACK_ILL_INPUT, "CPLAPACK", "CPLapackDense", MSGLS_BAD_NVECTOR);
    return(CPLAPACK_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */
  linit  = cpLapackDenseInit;
  lsetup = cpLapackDenseSetup;
  lsolve = cpLapackDenseSolve;
  lfree  = cpLapackDenseFree;

  /* Get memory for CPLapackMemRec */
  cplapack_mem = NULL;
  cplapack_mem = (CPLapackMem) malloc(sizeof(CPLapackMemRec));
  if (cplapack_mem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGLS_MEM_FAIL);
    return(CPLAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  djacE = NULL;
  djacI = NULL;
  J_data = NULL;

  last_flag = CPLAPACK_SUCCESS;
  lsetup_exists = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = LapackAllocDenseMat(N, N);
  if (M == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGLS_MEM_FAIL);
    free(cplapack_mem);
    return(CPLAPACK_MEM_FAIL);
  }
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(cplapack_mem);
    return(CPLAPACK_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = LapackAllocDenseMat(N, N);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDense", MSGLS_MEM_FAIL);
      LapackFreeMat(M);
      LapackFreeArray(pivots);
      free(cplapack_mem);
      return(CPLAPACK_MEM_FAIL);
    }
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cplapack_mem;

  return(CPLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cp_linit, cp_lsetup, cp_lsolve, and cp_lfree fields in (*cpode_mem)
 * to be cpLapackBandInit, cpLapackBandSetup, cpLapackBandSolve, 
 * and cpLapackBandFree, respectively.  It allocates memory for a 
 * structure of type CPLapackBandMemRec and sets the cp_lmem field in 
 * (*cpode_mem) to the address of this structure.  It sets lsetup_exists 
 * in (*cpode_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The CPLapackBand return value is CPLAPACK_SUCCESS = 0, 
 * CPLAPACK_MEM_FAIL = -1, or CPLAPACK_ILL_INPUT = -2.
 *
 * NOTE: The CPLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int CPLapackBand(void *cpode_mem, int N, int mupper, int mlower)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackBand", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    cpProcessError(cp_mem, CPLAPACK_ILL_INPUT, "CPLAPACK", "CPLapackBand", MSGLS_BAD_NVECTOR);
    return(CPLAPACK_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */  
  linit  = cpLapackBandInit;
  lsetup = cpLapackBandSetup;
  lsolve = cpLapackBandSolve;
  lfree  = cpLapackBandFree;
  
  /* Get memory for CPLapackMemRec */
  cplapack_mem = NULL;
  cplapack_mem = (CPLapackMem) malloc(sizeof(CPLapackMemRec));
  if (cplapack_mem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGLS_MEM_FAIL);
    return(CPLAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_BAND;

  /* Set default Jacobian routine and Jacobian data */
  bjacE = NULL;
  bjacI = NULL;
  J_data = NULL;

  last_flag = CPLAPACK_SUCCESS;
  lsetup_exists = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in cplapack_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    cpProcessError(cp_mem, CPLAPACK_ILL_INPUT, "CPLAPACK", "CPLapackBand", MSGLS_BAD_SIZES);
    return(CPLAPACK_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = LapackAllocBandMat(N, mu, ml, smu);
  if (M == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGLS_MEM_FAIL);
    free(cplapack_mem);
    return(CPLAPACK_MEM_FAIL);
  }  
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(cplapack_mem);
    return(CPLAPACK_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = LapackAllocBandMat(N, mu, ml, smu);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackBand", MSGLS_MEM_FAIL);
      LapackFreeMat(M);
      LapackFreeArray(pivots);
      free(cplapack_mem);
      return(CPLAPACK_MEM_FAIL);
    }
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cplapack_mem;

  return(CPLAPACK_SUCCESS);
}




/*
 * -----------------------------------------------------------------
 * Optional I/O functions for 
 * -----------------------------------------------------------------
 */

/*
 * CPLapackSetJacFn specifies the (dense or band) Jacobian function.
 */
int CPLapackSetJacFn(void *cpode_mem, void *jac, void *jac_data)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackSetJacFn", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackSetJacFn", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapack_mem = (CPLapackMem) lmem;

  switch (mtype) {
  case LAPACK_DENSE:
    if (ode_type == CP_EXPL) djacE = (CPLapackDenseJacExplFn) jac;
    else                     djacI = (CPLapackDenseJacImplFn) jac;
    break;
  case LAPACK_BAND:
    if (ode_type == CP_EXPL) bjacE = (CPLapackBandJacExplFn) jac;
    else                     bjacI = (CPLapackBandJacImplFn) jac;
    break;
  }
  J_data = jac_data;

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackGetWorkSpace returns the length of workspace allocated for the
 * CPLAPACK linear solver.
 */
int CPLapackGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackGetWorkSpace", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackGetWorkSpace", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapack_mem = (CPLapackMem) lmem;

  switch (mtype) {
  case LAPACK_DENSE:
    if (ode_type == CP_EXPL) *lenrwLS = 2*n*n;
    else                     *lenrwLS = n*n;
    *leniwLS = n;
    break;
  case LAPACK_BAND:
    break;
  }

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CPLapackGetNumJacEvals(void *cpode_mem, long int *njevals)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackGetNumJacEvals", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackGetNumJacEvals", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapack_mem = (CPLapackMem) lmem;

  *njevals = nje;

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackGetNumFctEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CPLapackGetNumFctEvals(void *cpode_mem, long int *nfevalsLS)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackGetNumFctEvals", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackGetNumFctEvals", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapack_mem = (CPLapackMem) lmem;

  *nfevalsLS = nfeDQ;

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackGetReturnFlagName returns the name associated with a CPLAPACK
 * return value.
 */
char *CPLapackGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPLAPACK_SUCCESS:
    sprintf(name,"CPLAPACK_SUCCESS");
    break;   
  case CPLAPACK_MEM_NULL:
    sprintf(name,"CPLAPACK_MEM_NULL");
    break;
  case CPLAPACK_LMEM_NULL:
    sprintf(name,"CPLAPACK_LMEM_NULL");
    break;
  case CPLAPACK_ILL_INPUT:
    sprintf(name,"CPLAPACK_ILL_INPUT");
    break;
  case CPLAPACK_MEM_FAIL:
    sprintf(name,"CPLAPACK_MEM_FAIL");
    break;
  case CPLAPACK_JACFUNC_UNRECVR:
    sprintf(name,"CPLAPACK_JACFUNC_UNRECVR");
    break;
  case CPLAPACK_JACFUNC_RECVR:
    sprintf(name,"CPLAPACK_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CPLapackGetLastFlag returns the last flag set in a CPLAPACK function.
 */
int CPLapackGetLastFlag(void *cpode_mem, int *flag)
{
  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackGetLastFlag", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackGetLastFlag", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapack_mem = (CPLapackMem) lmem;

  *flag = last_flag;

  return(CPLAPACK_SUCCESS);
}

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR PROJECTION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CPLapackDenseProj
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module for the projection.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPLapackDenseProj will
 *       first test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */
int CPLapackDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type)
{
  CPodeMem cp_mem;
  CPLapackProjMem cplapackP_mem;
  int ier, mone = -1;
  realtype tmp;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackDenseProj", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPLAPACK_ILL_INPUT, "CPLAPACK", "CPLapackDenseProj", MSGLS_BAD_NVECTOR);
    return(CPLAPACK_ILL_INPUT);
  }

  /* Check if fact_type has a legal value */
  if ( (fact_type != CPLAPACK_LU) && 
       (fact_type != CPLAPACK_QR) && 
       (fact_type != CPLAPACK_SC) &&
       (fact_type != CPLAPACK_QRP) ) {
    cpProcessError(cp_mem, CPLAPACK_ILL_INPUT, "CPLAPACK", "CPLapackDenseProj", MSGLS_BAD_FACT);
    return(CPLAPACK_ILL_INPUT);
  }

  if (lfreeP !=NULL) lfreeP(cp_mem);

  /* Set the five function fields in cp_mem */
  linitP  = cpLapackDenseProjInit;
  lsetupP = cpLapackDenseProjSetup;
  lsolveP = cpLapackDenseProjSolve;
  lmultP  = cpLapackDenseProjMult;
  lfreeP  = cpLapackDenseProjFree;

  /* Get memory for CPLapackProjMemRec */
  cplapackP_mem = NULL;
  cplapackP_mem = (CPLapackProjMem) malloc(sizeof(CPLapackProjMemRec));
  if (cplapackP_mem == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
    return(CPLAPACK_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  djacP = NULL;
  JP_data = NULL;

  lsetupP_exists = TRUE;

  /* Allocate memory */
  G = NULL;
  pivotsP = NULL;
  K = NULL;
  beta = NULL;
  wrk = NULL;
     
  /* Allocate memory for G */
  G = LapackAllocDenseMat(Ny, Nc);
  if (G == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
    free(cplapackP_mem);
    return(CPLAPACK_MEM_FAIL);
  }
  savedG = LapackAllocDenseMat(Ny, Nc);
  if (savedG == NULL) {
    cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
    LapackFreeMat(G);
    free(cplapackP_mem);
    return(CPLAPACK_MEM_FAIL);
  }

  /* Allocate additional work space, depending on factorization */
  switch(fact_type) {

  case CPLAPACK_LU:
    /* Allocate space for pivotsP and K */
    pivotsP = LapackAllocIntArray(Nc);
    if (pivotsP == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);
    }
    K = LapackAllocDenseMat(Ny-Nc, Ny-Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeArray(pivotsP);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);
    }
    break;

  case CPLAPACK_QR:
    /* Allocate space for beta */
    beta = LapackAllocRealArray(Nc);
    if (beta == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);      
    }
    /* Find optimal length of work array */
    dgeqrf_f77(&Ny, &Nc, G->data, &Ny, beta, &tmp, &mone, &ier);
    /* Allocate space for wrk */
    len_wrk = (int)tmp;
    wrk = LapackAllocRealArray(len_wrk);
    if (wrk == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeArray(beta);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = LapackAllocDenseMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
        LapackFreeArray(wrk);
        LapackFreeArray(beta);
        LapackFreeMat(savedG);
        LapackFreeMat(G);
        free(cplapackP_mem);
        return(CPLAPACK_MEM_FAIL);
      }
    }
    break;

  case CPLAPACK_QRP:
    /* Allocate space for pivotsP */
    pivotsP = LapackAllocIntArray(Nc);
    if (pivotsP == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);
    }
    /* Allocate space for beta */
    beta = LapackAllocRealArray(Nc);
    if (beta == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeArray(pivotsP);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);      
    }
    /* Find optimal length of work array */
    dgeqp3_f77(&Ny, &Nc, G->data, &Ny, pivotsP, beta, &tmp, &mone, &ier);
    /* Allocate space for wrk */
    len_wrk = (int)tmp;
    wrk = LapackAllocRealArray(len_wrk);
    if (wrk == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeArray(beta);
      LapackFreeArray(pivotsP);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = LapackAllocDenseMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
        LapackFreeArray(wrk);
        LapackFreeArray(beta);
        LapackFreeArray(pivotsP);
        LapackFreeMat(savedG);
        LapackFreeMat(G);
        free(cplapackP_mem);
        return(CPLAPACK_MEM_FAIL);
      }
    }
    break;

  case CPLAPACK_SC:
    /* Allocate space for K = G * D^(-1) * G^T */
    K = LapackAllocDenseMat(Nc, Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPLAPACK_MEM_FAIL, "CPLAPACK", "CPLapackDenseProj", MSGLS_MEM_FAIL);
      LapackFreeMat(savedG);
      LapackFreeMat(G);
      free(cplapackP_mem);
      return(CPLAPACK_MEM_FAIL);
    }

    break;

  }

  /* Copy inputs into memory */
  nc    = Nc;        /* number of constraints */
  ny    = Ny;        /* number of states      */
  ftype = fact_type; /* factorization type    */

  /* Attach linear solver memory to integrator memory */
  lmemP = cplapackP_mem;

  return(CPLAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional I/O functions for projection
 * -----------------------------------------------------------------
 */

/*
 * CPLapackProjSetJacFn specifies the constraint Jacobian function.
 */
int CPLapackProjSetJacFn(void *cpode_mem, void *jacP, void *jacP_data)
{
  CPodeMem cp_mem;
  CPLapackProjMem cplapackP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackProjSetJacFn", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackProjSetJacFn", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapackP_mem = (CPLapackProjMem) lmemP;

  djacP = (CPLapackDenseProjJacFn) jacP;
  JP_data = jacP_data;

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackProjGetNumJacEvals returns the number of constraint Jacobian evaluations
 */
int CPLapackProjGetNumJacEvals(void *cpode_mem, long int *njPevals)
{
  CPodeMem cp_mem;
  CPLapackProjMem cplapackP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackProjGetNumJacEvals", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackProjGetNumJacEvals", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapackP_mem = (CPLapackProjMem) lmemP;

  *njPevals = njeP;

  return(CPLAPACK_SUCCESS);
}

/*
 * CPLapackProjGetNumFctEvals returns the number of constraint function
 * evaluations for computing the DQ constraint Jacobian. 
 */
int CPLapackProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS)
{
  CPodeMem cp_mem;
  CPLapackProjMem cplapackP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPLAPACK_MEM_NULL, "CPLAPACK", "CPLapackProjGetNumFctEvals", MSGLS_CPMEM_NULL);
    return(CPLAPACK_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPLAPACK_LMEM_NULL, "CPLAPACK", "CPLapackProjGetNumFctEvals", MSGLS_LMEM_NULL);
    return(CPLAPACK_LMEM_NULL);
  }
  cplapackP_mem = (CPLapackProjMem) lmemP;

  *ncevalsLS = nceDQ;

  return(CPLAPACK_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cpLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int cpLapackDenseInit(CPodeMem cp_mem)
{
  CPLapackMem cplapack_mem;

  cplapack_mem = (CPLapackMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;
  
  if (ode_type == CP_EXPL && djacE == NULL) {
    djacE = cpLapackDenseDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && djacI == NULL) {
    djacI = cpLapackDenseDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPLAPACK_SUCCESS;
  return(0);
}

/*
 * cpLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the dense LU factorization routine.
 */
static int cpLapackDenseSetup(CPodeMem cp_mem, int convfail,
                              N_Vector yP, N_Vector ypP, N_Vector fctP,
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPLapackMem cplapack_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

  cplapack_mem = (CPLapackMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPL_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPL_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    if (jok) {

      /* If jok = TRUE, use saved copy of J */
      *jcurPtr = FALSE;
      dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);

    } else {

      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      retval = djacE(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
      } else if (retval < 0) {
        cpProcessError(cp_mem, CPLAPACK_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseSetup", MSGLS_JACFUNC_FAILED);
        last_flag = CPLAPACK_JACFUNC_UNRECVR;
        return(-1);
      } else if (retval > 0) {
        last_flag = CPLAPACK_JACFUNC_RECVR;
        return(1);
      }

    }

    /* Scale J by - gamma */
    fact = -gamma;
    dscal_f77(&(M->ldata), &fact, M->data, &one);

    /* Add identity to get M = I - gamma*J*/
    cpLapackDenseAddI(M);

    break;

  case CP_IMPL:

    /* Call Jacobian function */
    nje++;
    retval = djacI(n, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      break;
    } else if (retval < 0) {
      cpProcessError(cp_mem, CPLAPACK_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseSetup", MSGLS_JACFUNC_FAILED);
      last_flag = CPLAPACK_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPLAPACK_JACFUNC_RECVR;
      return(1);
    }
  
    break;

  }

  /* Do LU factorization of M */
  dgetrf_f77(&n, &n, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * cpLapackDenseAddI overwrites the dense matrix A with I+A
 */
static void cpLapackDenseAddI(LapackMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = A->cols[j];
    col_j[j] += ONE;
  }
}

/*
 * cpLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int cpLapackDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPLapackMem cplapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cplapack_mem = (CPLapackMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &n, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1); 
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }
  
  last_flag = CPLAPACK_SUCCESS;
  return(0);
}

/*
 * cpLapackDenseFree frees memory specific to the dense linear solver.
 */
static void cpLapackDenseFree(CPodeMem cp_mem)
{
  CPLapackMem  cplapack_mem;

  cplapack_mem = (CPLapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  if (ode_type == CP_EXPL) LapackFreeMat(savedJ);
  free(cplapack_mem); 
  cplapack_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * cpLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int cpLapackBandInit(CPodeMem cp_mem)
{
  CPLapackMem cplapack_mem;

  cplapack_mem = (CPLapackMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  if (ode_type == CP_EXPL && bjacE == NULL) {
    bjacE = cpLapackBandDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && bjacI == NULL) {
    bjacI = cpLapackBandDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPLAPACK_SUCCESS;
  return(0);
}

/*
 * cpLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the band LU factorization routine.
 */
static int cpLapackBandSetup(CPodeMem cp_mem, int convfail, 
                             N_Vector yP, N_Vector ypP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CPLapackMem cplapack_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;

  cplapack_mem = (CPLapackMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPL_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPL_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    if (jok) {
      
      /* If jok = TRUE, use saved copy of J */
      *jcurPtr = FALSE;
      dcopy_f77(&(savedJ->ldata), savedJ->data, &one, M->data, &one);
      
    } else {
      
      /* If jok = FALSE, call jac routine for new J value */
      nje++;
      nstlj = nst;
      *jcurPtr = TRUE;

      retval = bjacE(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      if (retval == 0) {
        dcopy_f77(&(M->ldata), M->data, &one, savedJ->data, &one);
      } else if (retval < 0) {
        cpProcessError(cp_mem, CPLAPACK_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackBandSetup", MSGLS_JACFUNC_FAILED);
        last_flag = CPLAPACK_JACFUNC_UNRECVR;
        return(-1);
      } else if (retval > 0) {
        last_flag = CPLAPACK_JACFUNC_RECVR;
        return(1);
      }

    }
  
    /* Scale J by - gamma */
    fact = -gamma;
    dscal_f77(&(M->ldata), &fact, M->data, &one);

    /* Add identity to get M = I - gamma*J*/
    cpLapackBandAddI(M);

    break;

  case CP_IMPL:

    /* Call Jacobian function */
    nje++;
    retval = bjacI(n, mu, ml, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      break;
    } else if (retval < 0) {
      cpProcessError(cp_mem, CPLAPACK_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackBandSetup", MSGLS_JACFUNC_FAILED);
      last_flag = CPLAPACK_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CPLAPACK_JACFUNC_RECVR;
      return(+1);
    }

    break;

  }

  /* Do LU factorization of M */
  dgbtrf_f77(&n, &n, &ml, &mu, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * cpLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int cpLapackBandSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPLapackMem cplapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  cplapack_mem = (CPLapackMem) lmem;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &n, &ml, &mu, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&n, &fact, bd, &one); 
  }

  last_flag = CPLAPACK_SUCCESS;
  return(0);
}

/*
 * cpLapackBandFree frees memory specific to the band linear solver.
 */
static void cpLapackBandFree(CPodeMem cp_mem)
{
  CPLapackMem  cplapack_mem;

  cplapack_mem = (CPLapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  if (ode_type == CP_EXPL) LapackFreeMat(savedJ);
  free(cplapack_mem); 
  cplapack_mem = NULL;
}

/*
 * cpLapackBandAddI overwrites the banded matrix A with I+A
 */
static void cpLapackBandAddI(LapackMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = (A->cols)[j] + A->storage_mu;
    col_j[0] += ONE;
  }
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR PROJECTION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cpLapackDenseProjInit does remaining initializations specific to
 * the dense linear solver.
 */
static int cpLapackDenseProjInit(CPodeMem cp_mem)
{
  CPLapackProjMem cplapackP_mem;

  cplapackP_mem = (CPLapackProjMem) lmemP;
  
  njeP   = 0;
  nceDQ  = 0;
  nstljP = 0;
  
  if (djacP == NULL) {
    djacP = cpLapackDenseProjDQJac;
    JP_data = cp_mem;
  }  

  return(0);
}

/*
 * cpLapackDenseProjSetup does the setup operations for the dense 
 * linear solver.
 * It calls the Jacobian evaluation routine and, depending on ftype,
 * it performs various factorizations.
 */
static int cpLapackDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1)
{
  int ier;
  CPLapackProjMem cplapackP_mem;
  realtype *col_i, rim1, ri;
  int i, j, k, nd, one = 1;
  int retval;

  realtype coef_1 = ONE, coef_0 = ZERO;

  cplapackP_mem = (CPLapackProjMem) lmemP;

  nd = ny-nc;

  /* Call Jacobian function (G will contain the Jacobian transposed) */
  retval = djacP(nc, ny, tn, y, cy, G, JP_data, c_tmp1, c_tmp2);
  njeP++;
  if (retval < 0) {
    cpProcessError(cp_mem, CPLAPACK_JACFUNC_UNRECVR, "CPLAPACK", "cpLapackDenseProjSetup", MSGLS_JACFUNC_FAILED);
    return(-1);
  } else if (retval > 0) {
    return(1);
  }

  /* Save Jacobian before factorization for possible use by lmultP */
  dcopy_f77(&(G->ldata), G->data, &one, savedG->data, &one);

  /* Factorize G, depending on ftype */
  switch (ftype) {

  case CPLAPACK_LU:

    /* 
     * LU factorization of G^T with partial pivoting
     *    P*G^T = | U1^T | * L^T
     *            | U2^T |
     * After factorization, P is encoded in pivotsP and
     * G^T is overwritten with U1 (nc by nc unit upper triangular), 
     * U2 ( nc by ny-nc rectangular), and L (nc by nc lower triangular).
     * Return ier if factorization failed. 
     */
    dgetrf_f77(&ny, &nc, G->data, &ny, pivotsP, &ier);
    if (ier != 0) return(ier);

    /* 
     * Build S = U1^{-1} * U2 (in place, S overwrites U2) 
     * For each row j of G, j = nc,...,ny-1, perform
     * a backward substitution (row version).
     * After this step, G^T contains U1, S, and L.
     */
    for (j=nc; j<ny; j++)
      dtrsv_f77("L", "T", "U", &nc, G->data, &ny, (G->data + j), &ny, 1, 1, 1);

    /*   
     * Build K = D1 + S^T * D2 * S 
     * S^T is stored in g_mat[nc...ny-1][0...nc]
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) {
      dsyrk_f77("L", "N", &nd, &nc, &coef_1, (G->data + nc), &ny, &coef_0, K->data, &nd, 1, 1);
      cpLapackDenseAddI(K);
    } else {
      cplLUcomputeKD(cp_mem, s_tmp1);
    }

    /*
     * Perform Cholesky decomposition of K: K = C*C^T
     * After factorization, the lower triangular part of K contains C.
     * Return ier if factorization failed. 
     */
    dpotrf_f77("L", &nd, K->data, &nd, &ier, 1);
    if (ier != 0) return(ier);

    break;

  case CPLAPACK_QR:

    /* 
     * QR factorization of G^T: G^T = Q*R
     * After factorization, the upper triangular part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     */
    dgeqrf_f77(&ny, &nc, G->data, &ny, beta, wrk, &len_wrk, &ier);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cplQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    break;

  case CPLAPACK_QRP:

    /* 
     * QR with pivoting factorization of G^T: G^T * P = Q * R.
     * After factorization, the upper triangular part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     * The pivots are stored in 'pivotsP'.
     */
    for (i=0; i<nc; i++) pivotsP[i] = 0;
    dgeqp3_f77(&ny, &nc, G->data, &ny, pivotsP, beta, wrk, &len_wrk, &ier);
    if (ier != 0) return(ier);

    /*
     * Determine the number of independent constraints.
     * After the QR factorization, the diagonal elements of R should 
     * be in decreasing order of their absolute values.
     */
    rim1 = ABS(G->data[0]);
    for (i=1, nr=1; i<nc; i++, nr++) {
      col_i = G->cols[i];
      ri = ABS(col_i[i]);
      if (ri < 100*uround) break;
      if (ri/rim1 < RPowerR(uround, THREE/FOUR)) break;
    }

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cplQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    break;

  case CPLAPACK_SC:

    /* 
     * Build K = G*D^(-1)*G^T
     * G^T is stored in g_mat[0...ny-1][0...nc]
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) {
      dsyrk_f77("L", "T", &nc, &ny, &coef_1, G->data, &ny, &coef_0, K->data, &nc, 1, 1);
    } else {
      cplSCcomputeKD(cp_mem, s_tmp1);
    }

    /* 
     * Perform Cholesky decomposition of K: K = C*C^T
     * After factorization, the lower triangular part of K contains C.
     * Return 1 if factorization failed. 
     */
    dpotrf_f77("L", &nc, K->data, &nc, &ier, 1);
    if (ier != 0) return(ier);

    break;

  }

  return(0);
}

/*
 * cpLapackDenseProjSolve handles the solve operation for the dense linear 
 * linear solver.
 */
static int cpLapackDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                                  N_Vector y, N_Vector cy,
                                  N_Vector c_tmp1, N_Vector s_tmp1)
{
  CPLapackProjMem cplapackP_mem;
  realtype *bd, *xd;
  realtype  *ewt_data, *d_data, *da_data, tmp;
  int nd, i, j, k, pk, ier, one = 1;
  realtype coef_1 = ONE, coef_0 = ZERO, coef_m1 = -ONE;

  cplapackP_mem = (CPLapackProjMem) lmemP;

  ewt_data = N_VGetArrayPointer(ewt);
  bd = N_VGetArrayPointer(b);
  xd = N_VGetArrayPointer(x);
  d_data = N_VGetArrayPointer(s_tmp1);

  nd = ny - nc;

  /* Solve the linear system, depending on ftype */
  switch (ftype) {

  case CPLAPACK_LU:

    /* Solve L*U1*alpha = bd
     *   (a) solve L*beta = bd using fwd. subst.
     *   (b) solve U1*alpha = beta using bckwd. subst
     * where L^T and U1^T are stored in G[0...nc-1][0...nc-1].
     * beta and then alpha overwrite bd.
     */
    dtrsv_f77("U", "T", "N", &nc, G->data, &ny, bd, &one, 1, 1, 1);
    dtrsv_f77("L", "T", "U", &nc, G->data, &ny, bd, &one, 1, 1, 1);

    /* Compute S^T * (D1 * alpha)
     * alpha is stored in bd.
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x2 = x[nc...ny-1].
     */
    if (pnorm == CP_PROJ_ERRNORM) {
      
      /* Load squared error weights into d */
      for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];
      /* Permute elements of d, based on pivotsP. 
       * Note that pivot information from dgetrf is 1-based.
       * Swap d[k] and d[pivotsP[k]]. */
      for (k=0; k<nc; k++) {
        pk = pivotsP[k];
        if (pk != k) {
          tmp = d_data[k];
          d_data[k]  = d_data[pk];
          d_data[pk] = tmp;
        }
      }
      /* Compute D1*alpha and store it into da_data */
      da_data = N_VGetArrayPointer(c_tmp1);
      for(k=0; k<nc; k++) da_data[k] = d_data[k] * bd[k];
      /* Compute S^T * D1 * alpha = S^T * da */
      dgemv_f77("N", &ny, &nc, &coef_1, (G->data + nc), &ny, da_data, &one, &coef_0, (xd + nc), &one, 1);
      
    } else {
      
      /* Compute S^T * alpha */
      dgemv_f77("N", &nd, &nc, &coef_1, (G->data + nc), &ny, bd, &one, &coef_0, (xd + nc), &one, 1);

    }

    /* Solve K*x2 = S^T*D1*alpha, using the Cholesky decomposition available in K.
     * S^T*D1*alpha is stored in x2 = x[nc...ny-1].
     */
    dpotrs_f77("L", &nd, &one, K->data, &nd, (xd + nc), &nd, &ier, 1);
    if (ier != 0) return(ier);

    /* Compute x1 = alpha - S*x2 
     * alpha is stored in bd.
     * x2 is stored in x[nc...ny-1].
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x1 = x[0...nc-1].
     */
    dcopy_f77(&nc, bd, &one, xd, &one);
    dgemv_f77("T", &nd, &nc, &coef_m1, (G->data + nc), &ny, (xd + nc), &one, &coef_1, xd, &one, 1);

    /* Compute P^T * x, where P is encoded into pivotsP.
     * Note that pivot information from dgetrf is 1-based.
     * Store result in x.
     */
    for (k=nc-1; k>=0; k--) {
      pk = pivotsP[k]-1;
      if(pk != k) {
        tmp = xd[k];
        xd[k] = xd[pk];
        xd[pk] = tmp;
      }
    }

    break;

  case CPLAPACK_QR:

    /* 
     * Solve R^T*alpha = bd using fwd. subst. (row version)
     * The upper triangular matrix R is stored in g_mat[0...nc-1][0...nc-1]
     * alpha overwrites bd.
     */
    dtrsv_f77("U", "T", "N", &nc, G->data, &ny, bd, &one, 1, 1, 1);

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) {
      dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    /* 
     * Compute x = Q1*alpha 
     *
     * Since we cannot really use the "thin" QR decomposition, we
     * first need to initialize xd = [alpha; 0].
     */
    for (k=0; k<nc; k++)  xd[k] = bd[k];
    for (k=nc; k<ny; k++) xd[k] = ZERO;
    dormqr_f77("L", "N", &ny, &one, &nc, G->data, &ny, beta, xd, &ny, wrk, &len_wrk, &ier, 1, 1);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;    


  case CPLAPACK_QRP:

    /* Compute P^T * b, where P is encoded into pivotsP.
     * If pivotsP[j] = k, then, for the factorization, the j-th column of G^T*P 
     * was the k-th column of G^T.
     * Therefore, to compute P^T*b, we must move the k-th element of b to the
     * j-th position, for j=1,2,... This is a forward permutation.
     * Note that pivot information from dgeqp3 is 1-based.
     * Store result in b.
     */
    for (i=1; i<=nc; i++) pivotsP[i-1] = -pivotsP[i-1];

    for (i=1; i<=nc; i++) {
      
      if (pivotsP[i-1] > 0) continue;

      j = i;
      pivotsP[j-1] = -pivotsP[j-1];
      pk = pivotsP[j-1];

      while (pivotsP[pk-1] < 0) {

        tmp = bd[j-1];
        bd[j-1] = bd[pk-1];
        bd[pk-1] = tmp;

        pivotsP[pk-1] = -pivotsP[pk-1];
        j = pk;
        pk = pivotsP[pk-1];
      }
      
    }

    /* 
     * Solve R11^T * alpha = P^T * bd using fwd. subst. (row version)
     * The upper triangular matrix R is stored in g_mat[0...nr-1][0...nr-1]
     * P^T * bd is available in bd.
     * We only consider the first nr components in P^T*bd.
     * alpha overwrites bd.
     */
    dtrsv_f77("U", "T", "N", &nr, G->data, &ny, bd, &one, 1, 1, 1);

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) {
      dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
      if (ier != 0) return(ier);
    }

    /* 
     * Compute x = Q1*alpha 
     *
     * Since we cannot really use the "thin" QR decomposition, we
     * first need to initialize xd = [alpha; 0].
     */
    for (k=0; k<nr; k++)  xd[k] = bd[k];
    for (k=nr; k<ny; k++) xd[k] = ZERO;
    dormqr_f77("L", "N", &ny, &one, &nc, G->data, &ny, beta, xd, &ny, wrk, &len_wrk, &ier, 1, 1);
    if (ier != 0) return(ier);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;    


  case CPLAPACK_SC:

    /* 
     * Solve K*xi = bd, using the Cholesky decomposition available in K.
     * xi overwrites bd.
     */
    dpotrs_f77("L", &nc, &one, K->data, &nc, bd, &nc, &ier, 1);
    if (ier != 0) return(ier);

    /* Compute x = G^T * xi
     * G^T is stored in g_mat[0...ny-1][0...nc-1]
     * xi is available in bd.
     */
    dgemv_f77("N", &ny, &nc, &coef_1, G->data, &ny, bd, &one, &coef_0, xd, &one, 1);

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) {
      for (i=0; i<ny; i++)
        xd[i] /= ewt_data[i]*ewt_data[i];
    }

    break;

  }

  return(0);
}

/*
 * cpDenseProjMult computes the Jacobian-vector product used a saved 
 * Jacobian copy.
 */
static void cpLapackDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx)
{
  CPLapackProjMem cplapackP_mem;
  realtype coef_1 = ONE, coef_0 = ZERO;
  realtype *xd, *Gxd;
  int one = 1;

  cplapackP_mem = (CPLapackProjMem) lmemP;

  xd = N_VGetArrayPointer(x);
  Gxd = N_VGetArrayPointer(Gx);

  dgemv_f77("T", &ny, &nc, &coef_1, savedG->data, &ny, xd, &one, &coef_0, Gxd, &one, 1);
}

/*
 * cpLapackDenseProjFree frees memory specific to the dense linear solver.
 */
static void cpLapackDenseProjFree(CPodeMem cp_mem)
{
  CPLapackProjMem cplapackP_mem;

  cplapackP_mem = (CPLapackProjMem) lmemP;
  
  LapackFreeMat(G);
  LapackFreeMat(savedG);
  switch (ftype) {
  case CPLAPACK_LU:
    LapackFreeArray(pivotsP);
    LapackFreeMat(K);
    break;
  case CPLAPACK_QR:
    LapackFreeArray(wrk);
    LapackFreeArray(beta);
    if (pnorm == CP_PROJ_ERRNORM) LapackFreeMat(K);
    break;
  case CPLAPACK_QRP:
    LapackFreeArray(wrk);
    LapackFreeArray(beta);
    LapackFreeArray(pivotsP);
    if (pnorm == CP_PROJ_ERRNORM) LapackFreeMat(K);
    break;
  case CPLAPACK_SC:
    LapackFreeMat(K);
    break;
  }

  free(cplapackP_mem); 
  cplapackP_mem = NULL;
}

/*
 * -----------------------------------------------------------------
 * Private functions for LU-, QR-, and SC-based projection
 * -----------------------------------------------------------------
 */

/*
 * Compute the lower triangle of K = D1 + S^T*D2*S,
 * D = diag(D1, D2) = P*W*P^T, W is a diagonal matrix
 * containing the squared error weights, and P is the 
 * permutation matrix encoded into pivotsP.
 * D1 has length nc and D2 has length (ny-nc).
 */
static void cplLUcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

static void cplQRcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

static void cplSCcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

/* 
 * =================================================================
 *  DENSE DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpLapackDenseDQJacExpl 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * LapackMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro LAPACK_DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
static int cpLapackDenseDQJacExpl(int N, realtype t,
                                  N_Vector y, N_Vector fy, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPLapackMem  cplapack_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cplapack_mem = (CPLapackMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(LAPACK_DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = fe(t, y, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    LAPACK_DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}


/*
 * -----------------------------------------------------------------
 * cpLapackDenseDQJacImpl 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian F_y' + gamma*F_y. It assumes that a dense matrix of type
 * LapackMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro LAPACK_DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
static int cpLapackDenseDQJacImpl(int N, realtype t, realtype gm,
                                  N_Vector y, N_Vector yp, N_Vector r, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data;
  N_Vector ftemp, jthCol;
  int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPLapackMem  cplapack_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cplapack_mem = (CPLapackMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y, and yp */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);
  yp_data  = N_VGetArrayPointer(yp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);

  /* Generate each column of the Jacobian M = F_y' + gamma * F_y
     as delta(F)/delta(y_j). */
  for (j = 0; j < N; j++) {

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(LAPACK_DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as h*yp_j. */
    inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewt_data[j] );
    if (h*ypj < ZERO) inc = -inc;
    inc = (yj + inc) - yj;

    /* Increment y_j and yp_j, call res, and break on error return. */
    y_data[j]  += gamma*inc;
    yp_data[j] += inc;
    retval = fi(t, y, yp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Generate the jth col of J(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, r, jthCol);

    LAPACK_DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

    /* Reset y_j, yp_j */     
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}

/* 
 * =================================================================
 *  BAND DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

static int cpLapackBandDQJacExpl(int N, int mupper, int mlower,
                                 realtype t, N_Vector y, N_Vector fy, 
                                 LapackMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int group, i, j, width, ngroups, i1, i2;
  int retval = 0;

  CPodeMem cp_mem;
  CPLapackMem  cplapack_mem;

  /* jac_dat points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cplapack_mem = (CPLapackMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */

    retval = fe(tn, ytemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = LAPACK_BAND_COL(Jac,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        LAPACK_BAND_COL_ELEM(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);
}

static int cpLapackBandDQJacImpl(int N, int mupper, int mlower,
                                 realtype t, realtype gm, 
                                 N_Vector y, N_Vector yp, N_Vector r,
                                 LapackMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  N_Vector ftemp, ytemp, yptemp;
  realtype inc, inc_inv, yj, ypj, srur, ewtj;
  realtype *y_data, *yp_data, *ewt_data;
  realtype *ytemp_data, *yptemp_data, *ftemp_data, *r_data, *col_j;
  int i, j, i1, i2, width, group, ngroups;
  int retval = 0;

  CPodeMem cp_mem;
  CPLapackMem cplapack_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cplapack_mem = (CPLapackMem) lmem;

  ftemp = tmp1; /* Rename work vector for use as the perturbed residual. */
  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */
  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

  /* Obtain pointers to the data for all vectors used.  */
  ewt_data    = N_VGetArrayPointer(ewt);
  r_data      = N_VGetArrayPointer(r);
  y_data      = N_VGetArrayPointer(y);
  yp_data     = N_VGetArrayPointer(yp);
  ftemp_data  = N_VGetArrayPointer(ftemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  /* Initialize ytemp and yptemp. */
  N_VScale(ONE, y, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */
  srur = RSqrt(uround);
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y[j] and yp[j] for j in this group. */
    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
           adjustments using ypj and ewtj if this is small, and a further
           adjustment to give it the same sign as h*ypj. */
        inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );

        if (h*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Increment yj and ypj. */
        ytemp_data[j]  += gamma*inc;
        yptemp_data[j] += inc;
    }

    /* Call ODE fct. with incremented arguments. */
    retval = fi(tn, ytemp, yptemp, ftemp, f_data);
    nfeDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */
    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */
      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = LAPACK_BAND_COL(Jac, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */
      inc = MAX( srur * MAX( ABS(yj), ABS(h*ypj) ) , ONE/ewtj );
      if (h*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      
      /* Load the difference quotient Jacobian elements for column j. */
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower,N-1);
      for (i=i1; i<=i2; i++) 
        LAPACK_BAND_COL_ELEM(col_j,i,j) = inc_inv*(ftemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);
}

/* 
 * =================================================================
 *  DENSE DQ JACOBIAN APPROXIMATIONS FOR PROJECTION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpLapackDenseProjDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation 
 * to the transpose of the Jacobian of c(t,y). It loads it into a 
 * dense matrix of type LapackMat stored column-wise with elements 
 * within each column contiguous. The address of the jth column of 
 * J is obtained via the macro LAPACK_DENSE_COL and this pointer is 
 * associated with an N_Vector using the N_VGetArrayPointer and 
 * N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian
 * transposed is done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
static int cpLapackDenseProjDQJac(int Nc, int Ny, realtype t,
                                  N_Vector y, N_Vector cy, 
                                  LapackMat Jac, void *jac_data,
                                  N_Vector c_tmp1, N_Vector c_tmp2)
{
  realtype inc, inc_inv, yj, srur;
  realtype *y_data, *ewt_data, *jthCol_data;
  N_Vector ctemp, jthCol;
  int i, j;
  int retval = 0;

  CPodeMem cp_mem;
  CPLapackProjMem cplapackP_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cplapackP_mem = (CPLapackProjMem) lmemP;

  /* Rename work vectors for readibility */
  ctemp  = c_tmp1; 
  jthCol = c_tmp2;

  /* Obtain pointers to the data for ewt and y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Obtain pointer to the data for jthCol */
  jthCol_data = N_VGetArrayPointer(jthCol);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);

  /* Generate each column of the Jacobian G = dc/dy as delta(c)/delta(y_j). */
  for (j = 0; j < Ny; j++) {

    /* Save the y_j values. */
    yj = y_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), 
       with an adjustment using ewt_j if this is small */
    inc = MAX( srur * ABS(yj) , ONE/ewt_data[j] );
    inc = (yj + inc) - yj;

    /* Increment y_j, call cfun, and break on error return. */
    y_data[j]  += inc;
    retval = cfun(t, y, ctemp, c_data);
    nceDQ++;
    if (retval != 0) break;

    /* Generate the jth col of G(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ctemp, -inc_inv, cy, jthCol);

    /* Copy the j-th column of G into the j-th row of Jac */
    for (i = 0; i < Nc ; i++) {
      LAPACK_DENSE_ELEM(Jac,j,i) = jthCol_data[i];
    }

    /* Reset y_j */     
    y_data[j] = yj;
  }

  return(retval);

}


