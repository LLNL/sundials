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
 * This is the implementation file for the CPDENSE linear solver.
 * -----------------------------------------------------------------
 */


/*
 * TODO:
 *
 *   cpdQRcomputeKD
 *   cpdSCcomputeKD
 *   cpdScaleXbyD
 */


/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_dense_impl.h"
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

/* CPDENSE linit, lsetup, lsolve, and lfree routines */
 
static int cpDenseInit(CPodeMem cp_mem);
static int cpDenseSetup(CPodeMem cp_mem, int convfail, 
                        N_Vector yP, N_Vector ypP, N_Vector fctP, 
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC);
static void cpDenseFree(CPodeMem cp_mem);

/* CPDENSE linitP, lsetupP, lsolveP, lmultP, and lfreeP routines */
static int cpDenseProjInit(CPodeMem cp_mem);
static int cpDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1);
static int cpDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                            N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector s_tmp1);
static void cpDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx);
static void cpDenseProjFree(CPodeMem cp_mem);

/* Private functions for LU, QR, and SC projection */
static void cpdLUcomputeKI(CPodeMem cp_mem);
static void cpdLUcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cpdQRcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cpdSCcomputeKI(CPodeMem cp_mem);
static void cpdSCcomputeKD(CPodeMem cp_mem, N_Vector d);
static void cpdScaleXbyD(CPodeMem cp_mem, realtype *alpha, realtype *d);

/* CPDENSE DQ integration Jacobian functions */
static int cpDenseDQJacExpl(long int N, realtype t,
                            N_Vector y, N_Vector fy, 
                            DenseMat Jac, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cpDenseDQJacImpl(long int N, realtype t, realtype gm,
                            N_Vector y, N_Vector yp, N_Vector r, 
                            DenseMat Jac, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* CPDENSE DQ projection Jacobian functions */
static int cpDenseProjDQJac(long int Nc, long int Ny, realtype t,
                            N_Vector y, N_Vector cy, 
                            DenseMat Jac, void *jac_data,
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
#define tempv          (cp_mem->cp_tempv)

#define linit          (cp_mem->cp_linit)
#define lsetup         (cp_mem->cp_lsetup)
#define lsolve         (cp_mem->cp_lsolve)
#define lfree          (cp_mem->cp_lfree)
#define lmem           (cp_mem->cp_lmem)
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

#define n              (cpdense_mem->d_n)
#define jacE           (cpdense_mem->d_jacE)
#define jacI           (cpdense_mem->d_jacI)
#define M              (cpdense_mem->d_M)
#define savedJ         (cpdense_mem->d_savedJ)
#define pivots         (cpdense_mem->d_pivots)
#define nstlj          (cpdense_mem->d_nstlj)
#define nje            (cpdense_mem->d_nje)
#define nfeD           (cpdense_mem->d_nfeD)
#define J_data         (cpdense_mem->d_J_data)
#define last_flag      (cpdense_mem->d_last_flag)

#define nc             (cpdenseP_mem->d_nc)
#define ny             (cpdenseP_mem->d_ny)
#define jacP           (cpdenseP_mem->d_jacP)
#define JP_data        (cpdenseP_mem->d_JP_data)
#define ftype          (cpdenseP_mem->d_ftype)
#define G              (cpdenseP_mem->d_G)
#define savedG         (cpdenseP_mem->d_savedG)
#define K              (cpdenseP_mem->d_K)
#define pivotsP        (cpdenseP_mem->d_pivotsP)
#define beta           (cpdenseP_mem->d_beta)
#define nstljP         (cpdenseP_mem->d_nstljP)
#define njeP           (cpdenseP_mem->d_njeP)
#define nceD           (cpdenseP_mem->d_nceD)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CPDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  CPDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cp_linit, cp_lsetup, cp_lsolve, cp_lfree fields in (*cpode_mem)
 * to be cpDenseInit, cpDenseSetup, cpDenseSolve, and cpDenseFree,
 * respectively.  It allocates memory for a structure of type
 * CPDenseMemRec and sets the cp_lmem field in (*cpode_mem) to the
 * address of this structure.  It sets lsetup_exists in (*cpode_mem) to
 * TRUE, and the d_jac field to the default cpDenseDQJac.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CPDense(void *cpode_mem, long int N)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDense", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDENSE_ILL_INPUT, "CPDENSE", "CPDense", MSGDS_BAD_NVECTOR);
    return(CPDENSE_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cp_mem);

  /* Set four main function fields in cp_mem */
  linit  = cpDenseInit;
  lsetup = cpDenseSetup;
  lsolve = cpDenseSolve;
  lfree  = cpDenseFree;

  /* Get memory for CPDenseMemRec */
  cpdense_mem = NULL;
  cpdense_mem = (CPDenseMem) malloc(sizeof(CPDenseMemRec));
  if (cpdense_mem == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDense", MSGDS_MEM_FAIL);
    return(CPDENSE_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  jacE = NULL;
  jacI = NULL;
  J_data = NULL;

  last_flag = CPDENSE_SUCCESS;
  lsetup_exists = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = DenseAllocMat(N, N);
  if (M == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDense", MSGDS_MEM_FAIL);
    free(cpdense_mem);
    return(CPDENSE_MEM_FAIL);
  }
  pivots = DenseAllocPiv(N);
  if (pivots == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDense", MSGDS_MEM_FAIL);
    DenseFreeMat(M);
    free(cpdense_mem);
    return(CPDENSE_MEM_FAIL);
  }
  if (ode_type == CP_EXPL) {
    savedJ = DenseAllocMat(N, N);
    if (savedJ == NULL) {
      cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDense", MSGDS_MEM_FAIL);
      DenseFreeMat(M);
      DenseFreePiv(pivots);
      free(cpdense_mem);
      return(CPDENSE_MEM_FAIL);
    }
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cpdense_mem;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseSetJacFn
 * -----------------------------------------------------------------
 */

int CPDenseSetJacFn(void *cpode_mem, void *djac, void *jac_data)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseSetJacFn", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseSetJacFn", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdense_mem = (CPDenseMem) lmem;

  if (ode_type == CP_EXPL) jacE = (CPDenseJacExplFn) djac;
  else                     jacI = (CPDenseJacImplFn) djac;
  J_data = jac_data;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetWorkSpace
 * -----------------------------------------------------------------
 */

int CPDenseGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseGetWorkSpace", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseGetWorkSpace", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdense_mem = (CPDenseMem) lmem;

  if (ode_type == CP_EXPL) *lenrwLS = 2*n*n;
  else                     *lenrwLS = n*n;
  *leniwLS = n;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetNumJacEvals
 * -----------------------------------------------------------------
 */

int CPDenseGetNumJacEvals(void *cpode_mem, long int *njevals)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseGetNumJacEvals", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseGetNumJacEvals", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdense_mem = (CPDenseMem) lmem;

  *njevals = nje;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetNumFctEvals
 * -----------------------------------------------------------------
 */

int CPDenseGetNumFctEvals(void *cpode_mem, long int *nfevalsLS)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseGetNumFctEvals", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseGetNumFctEvals", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdense_mem = (CPDenseMem) lmem;

  *nfevalsLS = nfeD;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *CPDenseGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CPDENSE_SUCCESS:
    sprintf(name,"CPDENSE_SUCCESS");
    break;   
  case CPDENSE_MEM_NULL:
    sprintf(name,"CPDENSE_MEM_NULL");
    break;
  case CPDENSE_LMEM_NULL:
    sprintf(name,"CPDENSE_LMEM_NULL");
    break;
  case CPDENSE_ILL_INPUT:
    sprintf(name,"CPDENSE_ILL_INPUT");
    break;
  case CPDENSE_MEM_FAIL:
    sprintf(name,"CPDENSE_MEM_FAIL");
    break;
  case CPDENSE_JACFUNC_UNRECVR:
    sprintf(name,"CPDENSE_JACFUNC_UNRECVR");
    break;
  case CPDENSE_JACFUNC_RECVR:
    sprintf(name,"CPDENSE_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetLastFlag
 * -----------------------------------------------------------------
 */

int CPDenseGetLastFlag(void *cpode_mem, int *flag)
{
  CPodeMem cp_mem;
  CPDenseMem cpdense_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseGetLastFlag", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmem == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseGetLastFlag", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdense_mem = (CPDenseMem) lmem;

  *flag = last_flag;

  return(CPDENSE_SUCCESS);
}

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR PROJECTION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * CPDenseProj
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module for the projection.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CPDenseProj(void *cpode_mem, long int Nc, long int Ny, int fact_type)
{
  CPodeMem cp_mem;
  CPDenseProjMem cpdenseP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseProj", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cpProcessError(cp_mem, CPDENSE_ILL_INPUT, "CPDENSE", "CPDenseProj", MSGDS_BAD_NVECTOR);
    return(CPDENSE_ILL_INPUT);
  }

  /* Check if fact_type has a legal value */
  if ( (fact_type != CPDENSE_LU) && (fact_type != CPDENSE_QR) && (fact_type != CPDENSE_SC) ) {
    cpProcessError(cp_mem, CPDENSE_ILL_INPUT, "CPDENSE", "CPDenseProj", MSGDS_BAD_FACT);
    return(CPDENSE_ILL_INPUT);
  }

  if (lfreeP !=NULL) lfreeP(cp_mem);

  /* Set the five function fields in cp_mem */
  linitP  = cpDenseProjInit;
  lsetupP = cpDenseProjSetup;
  lsolveP = cpDenseProjSolve;
  lmultP  = cpDenseProjMult;
  lfreeP  = cpDenseProjFree;

  /* Get memory for CPDenseProjMemRec */
  cpdenseP_mem = NULL;
  cpdenseP_mem = (CPDenseProjMem) malloc(sizeof(CPDenseProjMemRec));
  if (cpdenseP_mem == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
    return(CPDENSE_MEM_FAIL);
  }

  lsetupP_exists = TRUE;

  /* Initialize all internal pointers to NULL */
  G = NULL;
  K = NULL;
  pivotsP = NULL;
  beta = NULL;

  /* Allocate memory for G and other work space */
  G = DenseAllocMat(Ny, Nc);
  if (G == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
    free(cpdenseP_mem);
    return(CPDENSE_MEM_FAIL);
  }
  savedG = DenseAllocMat(Ny, Nc);
  if (savedG == NULL) {
    cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
    DenseFreeMat(G);
    free(cpdenseP_mem);
    return(CPDENSE_MEM_FAIL);
  }

  /* Allocate additional work space, depending on factorization */
  switch(fact_type) {

  case CPDENSE_LU:
    /* Allocate space for pivotsP and K */
    pivotsP = DenseAllocPiv(Nc);
    if (pivotsP == NULL) {
      cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
      DenseFreeMat(savedG);
      DenseFreeMat(G);
      free(cpdenseP_mem);
      return(CPDENSE_MEM_FAIL);
    }
    K = DenseAllocMat(Ny-Nc, Ny-Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
      DenseFreePiv(pivotsP);
      DenseFreeMat(savedG);
      DenseFreeMat(G);
      free(cpdenseP_mem);
      return(CPDENSE_MEM_FAIL);
    }
    break;

  case CPDENSE_QR:
    /* Allocate space for beta */
    beta = DenseAllocBeta(Nc);
    if (beta == NULL) {
      cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
      DenseFreeMat(savedG);
      DenseFreeMat(G);
      free(cpdenseP_mem);
      return(CPDENSE_MEM_FAIL);      
    }
    /* If projecting in WRMS norm, allocate space for K=Q^T*D^(-1)*Q */
    if (pnorm == CP_PROJ_ERRNORM) {
      K = DenseAllocMat(Nc, Nc);
      if (K == NULL) {
        cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
        DenseFreeBeta(beta);
      DenseFreeMat(savedG);
        DenseFreeMat(G);
        free(cpdenseP_mem);
        return(CPDENSE_MEM_FAIL);
      }
    }
    break;

  case CPDENSE_SC:
    /* Allocate space for K = G * D^(-1) * G^T */
    K = DenseAllocMat(Nc, Nc);
    if (K == NULL) {
      cpProcessError(cp_mem, CPDENSE_MEM_FAIL, "CPDENSE", "CPDenseProj", MSGDS_MEM_FAIL);
      DenseFreeMat(savedG);
      DenseFreeMat(G);
      free(cpdenseP_mem);
      return(CPDENSE_MEM_FAIL);
    }

    break;

  }

  /* Set default Jacobian routine and Jacobian data */
  jacP = NULL;
  JP_data = NULL;

  lsetupP_exists = TRUE;

  /* Copy inputs into memory */
  nc    = Nc;        /* number of constraints */
  ny    = Ny;        /* number of states      */
  ftype = fact_type; /* factorization type    */

  /* Attach linear solver memory to integrator memory */
  lmemP = cpdenseP_mem;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseProjSetJacFn
 * -----------------------------------------------------------------
 */

int CPDenseProjSetJacFn(void *cpode_mem, void *djacP, void *jacP_data)
{
  CPodeMem cp_mem;
  CPDenseProjMem cpdenseP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseProjSetJacFn", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseProjSetJacFn", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdenseP_mem = (CPDenseProjMem) lmemP;

  if (djacP != NULL) {
    jacP = (CPDenseProjJacFn) djacP;
    JP_data = jacP_data;
  }

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseProjGetNumJacEvals
 * -----------------------------------------------------------------
 */

int CPDenseProjGetNumJacEvals(void *cpode_mem, long int *njPevals)
{
  CPodeMem cp_mem;
  CPDenseProjMem cpdenseP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseProjGetNumJacEvals", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseProjGetNumJacEvals", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdenseP_mem = (CPDenseProjMem) lmemP;

  *njPevals = njeP;

  return(CPDENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CPDenseGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CPDenseProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS)
{
  CPodeMem cp_mem;
  CPDenseProjMem cpdenseP_mem;

  /* Return immediately if cpode_mem is NULL */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CPDENSE_MEM_NULL, "CPDENSE", "CPDenseProjGetNumFctEvals", MSGDS_CPMEM_NULL);
    return(CPDENSE_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  if (lmemP == NULL) {
    cpProcessError(cp_mem, CPDENSE_LMEM_NULL, "CPDENSE", "CPDenseProjGetNumFctEvals", MSGDS_LMEM_NULL);
    return(CPDENSE_LMEM_NULL);
  }
  cpdenseP_mem = (CPDenseProjMem) lmemP;

  *ncevalsLS = nceD;

  return(CPDENSE_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpDenseInit(CPodeMem cp_mem)
{
  CPDenseMem cpdense_mem;

  cpdense_mem = (CPDenseMem) lmem;
  
  nje   = 0;
  nfeD  = 0;
  nstlj = 0;
  
  if (ode_type == CP_EXPL && jacE == NULL) {
    jacE = cpDenseDQJacExpl;
    J_data = cp_mem;
  } 
  
  if (ode_type == CP_IMPL && jacI == NULL) {
    jacI = cpDenseDQJacImpl;
    J_data = cp_mem;
  }

  last_flag = CPDENSE_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy (for explicit ODE only). In any case, it constructs 
 * the Newton matrix M = I - gamma*J or M = F_y' - gamma*F_y, updates 
 * counters, and calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int cpDenseSetup(CPodeMem cp_mem, int convfail,
                        N_Vector yP, N_Vector ypP, N_Vector fctP,
                        booleantype *jcurPtr,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  CPDenseMem cpdense_mem;
  int retval;

  cpdense_mem = (CPDenseMem) lmem;

  switch (ode_type) {

  case CP_EXPL:

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma/gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CPD_MSBJ) ||
      ((convfail == CP_FAIL_BAD_J) && (dgamma < CPD_DGMAX)) ||
      (convfail == CP_FAIL_OTHER);
    jok = !jbad;
    
    /* Test if it is enough to use a saved Jacobian copy */
    if (jok) {
      *jcurPtr = FALSE;
      DenseCopy(savedJ, M);
    } else {
      nstlj = nst;
      *jcurPtr = TRUE;
      DenseZero(M);
      retval = jacE(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
      nje++;
      if (retval < 0) {
        cpProcessError(cp_mem, CPDENSE_JACFUNC_UNRECVR, "CPDENSE", "cpDenseSetup", MSGDS_JACFUNC_FAILED);
        last_flag = CPDENSE_JACFUNC_UNRECVR;
        return(-1);
      }
      if (retval > 0) {
        last_flag = CPDENSE_JACFUNC_RECVR;
        return(1);
      }
      DenseCopy(M, savedJ);
    }
  
    /* Scale and add I to get M = I - gamma*J */
    DenseScale(-gamma, M);
    DenseAddI(M);

    break;

  case CP_IMPL:

    /* Initialize Jacobian to 0 and call Jacobian function */
    DenseZero(M);
    retval = jacI(n, tn, gamma, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
    nje++;
    if (retval < 0) {
      cpProcessError(cp_mem, CPDENSE_JACFUNC_UNRECVR, "CPDENSE", "cpDenseSetup", MSGDS_JACFUNC_FAILED);
      last_flag = CPDENSE_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CPDENSE_JACFUNC_RECVR;
      return(1);
    }
  
    break;

  }

  /* Do LU factorization of M */
  ier = DenseGETRF(M, pivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cpDenseSolve(CPodeMem cp_mem, N_Vector b, N_Vector weight,
                        N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  CPDenseMem cpdense_mem;
  realtype *bd;

  cpdense_mem = (CPDenseMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, pivots, bd);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm_type == CP_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = CPDENSE_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cpDenseFree(CPodeMem cp_mem)
{
  CPDenseMem  cpdense_mem;

  cpdense_mem = (CPDenseMem) lmem;
  
  DenseFreeMat(M);
  DenseFreePiv(pivots);
  if (ode_type == CP_EXPL) DenseFreeMat(savedJ);
  free(cpdense_mem); 
  cpdense_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR PROJECTION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseProjInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cpDenseProjInit(CPodeMem cp_mem)
{
  CPDenseProjMem cpdenseP_mem;

  cpdenseP_mem = (CPDenseProjMem) lmemP;
  
  njeP   = 0;
  nceD   = 0;
  nstljP = 0;
  
  if (jacP == NULL) {
    jacP = cpDenseProjDQJac;
    JP_data = cp_mem;
  }  

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the Jacobian evaluation routine and, depending on ftype,
 * it performs various factorizations.
 * -----------------------------------------------------------------
 */

static int cpDenseProjSetup(CPodeMem cp_mem, N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector c_tmp2, N_Vector s_tmp1)
{
  long int ier;
  CPDenseProjMem cpdenseP_mem;
  realtype **g_mat, *col_i, *s_tmpd;
  long int i, j, k, one = 1;
  int retval;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  g_mat = G->data;

  /* 
   * Initialize Jacobian matrix to 0 and call Jacobian function 
   * G will contain the Jacobian transposed. 
   */
  DenseZero(G);
  retval = jacP(nc, ny, tn, y, cy, G, JP_data, c_tmp1, c_tmp2);
  njeP++;
  if (retval < 0) {
    cpProcessError(cp_mem, CPDENSE_JACFUNC_UNRECVR, "CPDENSE", "cpDenseProjSetup", MSGDS_JACFUNC_FAILED);
    return(-1);
  } else if (retval > 0) {
    return(1);
  }

  /* Save Jacobian before factorization for possible use by lmultP */
  DenseCopy(G, savedG);

  /* Factorize G, depending on ftype */
  switch (ftype) {

  case CPDENSE_LU:

    /* 
     * LU factorization of G^T
     *      
     *    P*G^T =  | U1^T | * L^T     
     *             | U2^T |
     *
     * After factorization, P is encoded in pivotsP and
     * G^T is overwritten with U1 (nc by nc unit upper triangular), 
     * U2 ( nc by ny-nc rectangular), and L (nc by nc lower triangular).
     *
     * Return 1 if factorization failed. 
     */
    ier = DenseGETRF(G, pivotsP); 
    if (ier > 0) return(1);

    /* 
     * Build S = U1^{-1} * U2 (in place, S overwrites U2) 
     * For each row j of G, j = nc,...,ny-1, perform
     * a backward substitution (row version).
     *
     * After this step, G^T contains U1, S, and L.
     */
    for (j=nc; j<ny; j++) {
      for (i=nc-2; i>=0; i--) {
        col_i = g_mat[i];
        for (k=i+1; k<nc; k++) g_mat[i][j] -= col_i[k]*g_mat[k][j];
      }      
    }

    /* 
     * Build K = D1 + S^T * D2 * S 
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) cpdLUcomputeKI(cp_mem);
    else                         cpdLUcomputeKD(cp_mem, s_tmp1);

    /* 
     * Perform Cholesky decomposition of K (in place, gaxpy version)
     * 
     *     K = C*C^T
     *
     * After factorization, the lower triangular part of K contains C.
     *
     * Return 1 if factorization failed. 
     */
    ier = DensePOTRF(K);
    if (ier > 0) return(1);

    break;

  case CPDENSE_QR:

    /* 
     * Thin QR factorization of G^T
     *
     *   G^T = Q * R
     *
     * After factorization, the upper trianguler part of G^T 
     * contains the matrix R. The lower trapezoidal part of
     * G^T, together with the array beta, encodes the orthonormal
     * columns of Q as elementary reflectors.
     */

    /* Use s_tmp1 as workspace */
    s_tmpd = N_VGetArrayPointer(s_tmp1);
    ier = DenseGEQRF(G, beta, s_tmpd);

    /* If projecting in WRMS norm */
    if (pnorm == CP_PROJ_ERRNORM) {
      /* Build K = Q^T * D^(-1) * Q */
      cpdQRcomputeKD(cp_mem, s_tmp1);
      /* Perform Cholesky decomposition of K */
      ier = DensePOTRF(K);
      if (ier > 0) return(1);
    }

    break;

  case CPDENSE_SC:

    /* 
     * Build K = G*D^(-1)*G^T
     * Compute and store only the lower triangular part of K.
     */
    if (pnorm == CP_PROJ_L2NORM) cpdSCcomputeKI(cp_mem);
    else                         cpdSCcomputeKD(cp_mem, s_tmp1);

    /* 
     * Perform Cholesky decomposition of K (in place, gaxpy version)
     * 
     *     K = C*C^T
     *
     * After factorization, the lower triangular part of K contains C.
     *
     * Return 1 if factorization failed. 
     */
    ier = DensePOTRF(K);
    if (ier > 0) return(1);

    break;

  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear 
 * solver. The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cpDenseProjSolve(CPodeMem cp_mem, N_Vector b, N_Vector x,
                            N_Vector y, N_Vector cy,
                            N_Vector c_tmp1, N_Vector s_tmp1)
{
  CPDenseProjMem cpdenseP_mem;
  realtype **g_mat, *bd, *xd, *col_i, *s_tmpd;
  realtype  *ewt_data, *d_data, *da_data, tmp;
  long int nd, i, k, pk;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  g_mat = G->data;
  bd = N_VGetArrayPointer(b);
  xd = N_VGetArrayPointer(x);
  d_data = N_VGetArrayPointer(s_tmp1);

  nd = ny - nc;

  /* Solve the linear system, depending on ftype */
  switch (ftype) {

  case CPDENSE_LU:

    /* 
     * Solve L*U1*alpha = bd
     *   (a) solve L*beta = bd using fwd. subst. (row version)
     *   (b) solve U1*alpha = beta using bckwd. subst (row version) 
     * where L^T and U1^T are stored in G[0...nc-1][0...nc-1].
     * beta and then alpha overwrite bd.
     */
    bd[0] /= g_mat[0][0];
    for (i=1; i<nc; i++) {
      col_i = g_mat[i];
      for (k=0; k<i; k++) bd[i] -= col_i[k]*bd[k];
      bd[i] /= col_i[i];
    }
    for (i=nc-2; i>=0; i--) {
      col_i = g_mat[i];
      for (k=i+1; k<nc; k++) bd[i] -= col_i[k]*bd[k];
    }  

    /* 
     * Compute S^T * (D1 * alpha)
     * alpha is stored in bd.
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x2 = x[nc...ny-1].
     */
    if (pnorm == CP_PROJ_ERRNORM) {

      ewt_data = N_VGetArrayPointer(ewt);
      /* Load squared error weights into d */
      for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];
      /* Permute elements of d, based on pivotsP. Swap d[k] and d[pivotsP[k]]. */
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
      for(i=0; i<nd; i++) {
        xd[nc+i] = ZERO;
        for(k=0; k<nc; k++) xd[nc+i] += g_mat[k][nc+i]*da_data[k];
      }

    } else {

      /* Compute S^T * alpha */
      for(i=0; i<nd; i++) {
        xd[nc+i] = ZERO;
        for(k=0; k<nc; k++) xd[nc+i] += g_mat[k][nc+i]*bd[k];
      }

    }

    /* 
     * Solve K*x2 = S^T*D1*alpha, using the Cholesky decomposition available in K.
     * S^T*D1*alpha is stored in x2 = x[nc...ny-1].
     */
    DensePOTRS(K, &xd[nc]);

    /* 
     * Compute x1 = alpha - S*x2 
     * alpha is stored in bd.
     * x2 is stored in x[nc...ny-1].
     * S^T is stored in g_mat[nc...ny-1][0...nc-1].
     * Store result in x1 = x[0...nc-1].
     */
    for (i=0; i<nc; i++) {
      xd[i] = bd[i];
      col_i = g_mat[i];
      for (k=nc; k<ny; k++) xd[i] -= col_i[k]*xd[k];
    }

    /* 
     * Compute P^T * x, where P is encoded into pivotsP.
     * Store result in x.
     */
    for (k=nc-1; k>=0; k--) {
      pk = pivotsP[k];
      if(pk != k) {
        tmp = xd[k];
        xd[k] = xd[pk];
        xd[pk] = tmp;
      }
    }


    break;

  case CPDENSE_QR:

    /* 
     * Solve R^T*alpha = bd using fwd. subst. (row version)
     * alpha overwrites bd.
     */
    bd[0] /= g_mat[0][0];
    for (i=1; i<nc; i++) {
      col_i = g_mat[i];
      for (k=0; k<i; k++) bd[i] -= bd[k]*col_i[k];
      bd[i] /= col_i[i];
    }

    /* If projecting in WRMS norm, solve K*beta = alpha */
    if (pnorm == CP_PROJ_ERRNORM) DensePOTRS(K, bd);

    /* Compute x = Q*alpha */
    s_tmpd = N_VGetArrayPointer(s_tmp1);
    DenseORMQR(G, beta, bd, xd, s_tmpd); 

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) cpdScaleXbyD(cp_mem, xd, d_data);

    break;

  case CPDENSE_SC:

    /* 
     * Solve K*xi = bd, using the Cholesky decomposition available in K.
     * xi overwrites bd.
     */
    DensePOTRS(K, bd);

    /* Compute x = G^T * xi
     * G^T is stored in g_mat[0...ny-1][0...nc-1]
     */
    for(i=0; i<ny; i++) {
      xd[i] = ZERO;
      for(k=0; k<nc; k++) xd[i] += g_mat[k][i]*bd[k];
    }

    /* If projecting in WRMS norm, scale x by D^(-1) */
    if (pnorm == CP_PROJ_ERRNORM) cpdScaleXbyD(cp_mem, xd, d_data);

    break;

  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjMult
 * -----------------------------------------------------------------
 * This routine computes the Jacobian-vector product used a saved 
 * Jacobian copy.
 * -----------------------------------------------------------------
 */

static void cpDenseProjMult(CPodeMem cp_mem, N_Vector x, N_Vector Gx)
{
  CPDenseProjMem cpdenseP_mem;
  realtype **g_mat, *col_j, *x_d, *Gx_d;
  long int j, k;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  g_mat = savedG->data;
  x_d   = N_VGetArrayPointer(x);
  Gx_d  = N_VGetArrayPointer(Gx);

  for (j=0; j<nc; j++) {
    col_j = g_mat[j];
    Gx_d[j] = ZERO;
    for (k=0; k<ny; k++) Gx_d[j] += col_j[k]*x_d[k];
  }
}

/*
 * -----------------------------------------------------------------
 * cpDenseProjFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cpDenseProjFree(CPodeMem cp_mem)
{
  CPDenseProjMem  cpdenseP_mem;

  cpdenseP_mem = (CPDenseProjMem) lmemP;
  
  DenseFreeMat(G);
  DenseFreeMat(savedG);
  switch (ftype) {
  case CPDENSE_LU:
    DenseFreePiv(pivotsP);
    DenseFreeMat(K);
    break;
  case CPDENSE_QR:
    DenseFreeBeta(beta);
    if (pnorm == CP_PROJ_ERRNORM) DenseFreeMat(K);
    break;
  case CPDENSE_SC:
    DenseFreeMat(K);
    break;
  }

  free(cpdenseP_mem); 
  cpdenseP_mem = NULL;
}

/*
 * -----------------------------------------------------------------
 * Private functions for LU-, QR-, and SC-based projection
 * -----------------------------------------------------------------
 */

/*
 * Compute the lower triangle of K = I + S^T*S
 */
static void cpdLUcomputeKI(CPodeMem cp_mem)
{
  CPDenseProjMem cpdenseP_mem;
  realtype **g_mat, **k_mat, *k_col_j;
  long int nd, i, j, k;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  g_mat = G->data;
  k_mat = K->data;

  nd = ny-nc;

  /* Load K column by column */
  for (j=0; j<nd; j++) {

    k_col_j = k_mat[j];
    
    for (i=j; i<nd; i++) {
      k_col_j[i] = ZERO;
      for (k=0; k<nc; k++) k_col_j[i] += g_mat[k][nc+i]*g_mat[k][nc+j];
    }
    
    k_col_j[j] += ONE;

  }
  
}

/*
 * Compute the lower triangle of K = D1 + S^T*D2*S,
 * D = diag(D1, D2) = P*W*P^T, W is a diagonal matrix
 * containing the squared error weights, and P is the 
 * permutation matrix encoded into pivotsP.
 * D1 has length nc and D2 has length (ny-nc).
 */
static void cpdLUcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  CPDenseProjMem cpdenseP_mem;
  realtype *d_data, *ewt_data, tmp;
  realtype **g_mat, **k_mat, *k_col_j;
  long int nd, i, j, k, pk;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  ewt_data = N_VGetArrayPointer(ewt);
  d_data   = N_VGetArrayPointer(d);

  g_mat = G->data;
  k_mat = K->data;

  nd = ny-nc;

  /* Load squared error weights into d */
  for (k=0; k<ny; k++) d_data[k] = ewt_data[k] * ewt_data[k];

  /* Permute elements of d, based on pivotsP. Swap d[k] and d[pivotsP[k]]. */
  for (k=0; k<nc; k++) {
    pk = pivotsP[k];
    if (pk != k) {
      tmp = d_data[k];
      d_data[k]  = d_data[pk];
      d_data[pk] = tmp;
    }
  }

  /* load K column by column */
  for (j=0; j<nd; i++) {

    k_col_j = k_mat[j];
 
    for (i=j; i<nd; i++) {
      k_col_j[i] = ZERO;
      for(k=0; k<nc; k++) 
        k_col_j[i] += g_mat[k][nc+i] * d_data[k]*d_data[k] * g_mat[k][nc+j];
    }

    k_col_j[j] += d_data[j]*d_data[j];

  }

}

static void cpdQRcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}


/*
 * Compute the lower triangle of K = G * G^T
 * G^T is available in g_mat[0...ny][0...nc].
 */
static void cpdSCcomputeKI(CPodeMem cp_mem)
{
  CPDenseProjMem cpdenseP_mem;
  realtype **g_mat, **k_mat, *k_col_j;
  long int nd, i, j, k;

  cpdenseP_mem = (CPDenseProjMem) lmemP;

  g_mat = G->data;
  k_mat = K->data;

  /* Load K column by column */
  for (j=0; j<nc; j++) {
    k_col_j = k_mat[j];
    for (i=0; i<nc; i++) {
      k_col_j[i] = ZERO;
      for (k=0; k<ny; k++) k_col_j[i] += g_mat[i][k]*g_mat[j][k];
    }    
  }

}

static void cpdSCcomputeKD(CPodeMem cp_mem, N_Vector d)
{
  /* RADU:: implement this ... */
}

static void cpdScaleXbyD(CPodeMem cp_mem, realtype *alpha, realtype *d)
{
  /* RADU:: implement this ... */
}



/* 
 * =================================================================
 *  DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseDQJacExpl 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * DenseMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */
 
static int cpDenseDQJacExpl(long int N, realtype t,
                            N_Vector y, N_Vector fy, 
                            DenseMat Jac, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDenseMem  cpdense_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdense_mem = (CPDenseMem) lmem;

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

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = fe(t, y, ftemp, f_data);
    nfeD++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}


/*
 * -----------------------------------------------------------------
 * cpDenseDQJacImpl 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian F_y' + gamma*F_y. It assumes that a dense matrix of type
 * DenseMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */
 
static int cpDenseDQJacImpl(long int N, realtype t, realtype gm,
                            N_Vector y, N_Vector yp, N_Vector r, 
                            DenseMat Jac, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDenseMem  cpdense_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdense_mem = (CPDenseMem) lmem;

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
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);
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
    nfeD++;
    if (retval != 0) break;

    /* Generate the jth col of J(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, r, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

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
 *  DQ JACOBIAN APPROXIMATIONS FOR PROJECTION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cpDenseProjDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation 
 * to the transpose of the Jacobian of c(t,y). It loads it into a 
 * dense matrix of type DenseMat stored column-wise with elements 
 * within each column contiguous. The address of the jth column of 
 * J is obtained via the macro DENSE_COL and this pointer is 
 * associated with an N_Vector using the N_VGetArrayPointer and 
 * N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian
 * transposed is done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */
 
static int cpDenseProjDQJac(long int Nc, long int Ny, realtype t,
                            N_Vector y, N_Vector cy, 
                            DenseMat Jac, void *jac_data,
                            N_Vector c_tmp1, N_Vector c_tmp2)
{
  realtype inc, inc_inv, yj, srur;
  realtype *y_data, *ewt_data, *jthCol_data;
  N_Vector ctemp, jthCol;
  long int i, j;
  int retval = 0;

  CPodeMem cp_mem;
  CPDenseProjMem cpdenseP_mem;

  /* jac_data points to cpode_mem */
  cp_mem = (CPodeMem) jac_data;
  cpdenseP_mem = (CPDenseProjMem) lmemP;

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
    nceD++;
    if (retval != 0) break;

    /* Generate the jth col of G(tn,y) */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ctemp, -inc_inv, cy, jthCol);

    /* Copy the j-th column of G into the j-th row of Jac */
    for (i = 0; i < Nc ; i++) {
      DENSE_ELEM(Jac,j,i) = jthCol_data[i];
    }

    /* Reset y_j */     
    y_data[j] = yj;
  }

  return(retval);

}

