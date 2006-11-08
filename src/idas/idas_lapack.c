/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:28 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a IDAS dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "idas_lapack_impl.h"
#include "idas_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* IDALAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int idaLapackDenseInit(IDAMem IDA_mem);
static int idaLapackDenseSetup(IDAMem IDA_mem,
                               N_Vector yP, N_Vector ypP, N_Vector fctP, 
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int idaLapackDenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                               N_Vector yC, N_Vector ypC, N_Vector fctC);
static int idaLapackDenseFree(IDAMem IDA_mem);

/* IDALAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int idaLapackBandInit(IDAMem IDA_mem);
static int idaLapackBandSetup(IDAMem IDA_mem,
                              N_Vector yP, N_Vector ypP, N_Vector fctP, 
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int idaLapackBandSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC);
static int idaLapackBandFree(IDAMem IDA_mem);

/* IDALAPACK DENSE and BAND DQ integration Jacobian functions */
static int idaLapackDenseDQJac(int N, realtype tt, realtype c_j,
                               N_Vector yy, N_Vector yp, N_Vector rr, 
                               LapackMat Jac, void *jac_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int idaLapackBandDQJac(int N, int mupper, int mlower,
                              realtype tt, realtype c_j, 
                              N_Vector yy, N_Vector yp, N_Vector rr,
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define res            (IDA_mem->ida_res)
#define rdata          (IDA_mem->ida_rdata)
#define uround         (IDA_mem->ida_uround)
#define nst            (IDA_mem->ida_nst)
#define tn             (IDA_mem->ida_tn)
#define hh             (IDA_mem->ida_hh)
#define cj             (IDA_mem->ida_cj)
#define cjratio        (IDA_mem->ida_cjratio)
#define ewt            (IDA_mem->ida_ewt)
#define constraints    (IDA_mem->ida_constraints)

#define linit          (IDA_mem->ida_linit)
#define lsetup         (IDA_mem->ida_lsetup)
#define lsolve         (IDA_mem->ida_lsolve)
#define lfree          (IDA_mem->ida_lfree)
#define lperf          (IDA_mem->ida_lperf)
#define lmem           (IDA_mem->ida_lmem)
#define tempv          (IDA_mem->ida_tempv1)
#define setupNonNull   (IDA_mem->ida_setupNonNull)

#define mtype          (idalapack_mem->l_mtype)
#define n              (idalapack_mem->l_n)
#define ml             (idalapack_mem->b_ml)
#define mu             (idalapack_mem->b_mu)
#define smu            (idalapack_mem->b_smu)
#define djac           (idalapack_mem->d_jac)
#define bjac           (idalapack_mem->b_jac)
#define M              (idalapack_mem->l_M)
#define pivots         (idalapack_mem->l_pivots)
#define nje            (idalapack_mem->l_nje)
#define nreDQ          (idalapack_mem->l_nreDQ)
#define J_data         (idalapack_mem->l_J_data)
#define last_flag      (idalapack_mem->l_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * IDALapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  IDALapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the ida_linit, ida_lsetup, ida_lsolve, ida_lfree fields in (*ida_mem)
 * to be idaLapackDenseInit, idaLapackDenseSetup, idaLapackDenseSolve, 
 * and idaLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type IDALapackMemRec and sets the ida_lmem field in 
 * (*ida_mem) to the address of this structure.  It sets setupNonNull 
 * in (*ida_mem) to TRUE, and the d_jac field to the default 
 * idaLapackDenseDQJac. Finally, it allocates memory for M, pivots.
 *
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDALapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int IDALapackDense(void *ida_mem, int N)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackDense", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_ILL_INPUT, "IDALAPACK", "IDALapackDense", MSGLS_BAD_NVECTOR);
    return(IDALAPACK_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(IDA_mem);

  /* Set four main function fields in IDA_mem */
  linit  = idaLapackDenseInit;
  lsetup = idaLapackDenseSetup;
  lsolve = idaLapackDenseSolve;
  lperf  = NULL;
  lfree  = idaLapackDenseFree;

  /* Get memory for IDALapackMemRec */
  idalapack_mem = NULL;
  idalapack_mem = (IDALapackMem) malloc(sizeof(IDALapackMemRec));
  if (idalapack_mem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackDense", MSGLS_MEM_FAIL);
    return(IDALAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  djac = NULL;
  J_data = NULL;

  last_flag = IDALAPACK_SUCCESS;
  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, pivot array */
  M = NULL;
  pivots = NULL;

  M = LapackAllocDenseMat(N, N);
  if (M == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackDense", MSGLS_MEM_FAIL);
    free(idalapack_mem);
    return(IDALAPACK_MEM_FAIL);
  }
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackDense", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(idalapack_mem);
    return(IDALAPACK_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = idalapack_mem;

  return(IDALAPACK_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDALapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * ida_linit, ida_lsetup, ida_lsolve, and ida_lfree fields in (*ida_mem)
 * to be idaLapackBandInit, idaLapackBandSetup, idaLapackBandSolve, 
 * and idaLapackBandFree, respectively.  It allocates memory for a 
 * structure of type IDALapackBandMemRec and sets the ida_lmem field in 
 * (*ida_mem) to the address of this structure.  It sets setupNonNull 
 * in (*ida_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M and pivots.
 * The IDALapackBand return value is IDALAPACK_SUCCESS = 0, 
 * IDALAPACK_MEM_FAIL = -1, or IDALAPACK_ILL_INPUT = -2.
 *
 * NOTE: The IDALAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDALapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int IDALapackBand(void *ida_mem, int N, int mupper, int mlower)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackBand", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_ILL_INPUT, "IDALAPACK", "IDALapackBand", MSGLS_BAD_NVECTOR);
    return(IDALAPACK_ILL_INPUT);
  }

  if (lfree != NULL) lfree(IDA_mem);

  /* Set four main function fields in IDA_mem */  
  linit  = idaLapackBandInit;
  lsetup = idaLapackBandSetup;
  lsolve = idaLapackBandSolve;
  lperf  = NULL;
  lfree  = idaLapackBandFree;
  
  /* Get memory for IDALapackMemRec */
  idalapack_mem = NULL;
  idalapack_mem = (IDALapackMem) malloc(sizeof(IDALapackMemRec));
  if (idalapack_mem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackBand", MSGLS_MEM_FAIL);
    return(IDALAPACK_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = LAPACK_BAND;

  /* Set default Jacobian routine and Jacobian data */
  bjac = NULL;
  J_data = NULL;

  last_flag = IDALAPACK_SUCCESS;
  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in idalapack_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    IDAProcessError(IDA_mem, IDALAPACK_ILL_INPUT, "IDALAPACK", "IDALapackBand", MSGLS_BAD_SIZES);
    return(IDALAPACK_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M and pivot arrays */
  M = NULL;
  pivots = NULL;

  M = LapackAllocBandMat(N, mu, ml, smu);
  if (M == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackBand", MSGLS_MEM_FAIL);
    free(idalapack_mem);
    return(IDALAPACK_MEM_FAIL);
  }  
  pivots = LapackAllocIntArray(N);
  if (pivots == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_MEM_FAIL, "IDALAPACK", "IDALapackBand", MSGLS_MEM_FAIL);
    LapackFreeMat(M);
    free(idalapack_mem);
    return(IDALAPACK_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = idalapack_mem;

  return(IDALAPACK_SUCCESS);
}




/*
 * -----------------------------------------------------------------
 * Optional I/O functions for 
 * -----------------------------------------------------------------
 */

/*
 * IDALapackSetJacFn specifies the (dense or band) Jacobian function.
 */
int IDALapackSetJacFn(void *ida_mem, void *jac, void *jac_data)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackSetJacFn", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_LMEM_NULL, "IDALAPACK", "IDALapackSetJacFn", MSGLS_LMEM_NULL);
    return(IDALAPACK_LMEM_NULL);
  }
  idalapack_mem = (IDALapackMem) lmem;

  if (mtype == LAPACK_DENSE) 
    djac = (IDALapackDenseJacFn) jac;
  else if (mtype == LAPACK_BAND)
    bjac = (IDALapackBandJacFn) jac;

  J_data = jac_data;

  return(IDALAPACK_SUCCESS);
}

/*
 * IDALapackGetWorkSpace returns the length of workspace allocated for the
 * IDALAPACK linear solver.
 */
int IDALapackGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackGetWorkSpace", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_LMEM_NULL, "IDALAPACK", "IDALapackGetWorkSpace", MSGLS_LMEM_NULL);
    return(IDALAPACK_LMEM_NULL);
  }
  idalapack_mem = (IDALapackMem) lmem;

  if (mtype == LAPACK_DENSE) {
    *lenrwLS = n*n;
    *leniwLS = n;
  } else if (mtype == LAPACK_BAND) {

  }
    
  return(IDALAPACK_SUCCESS);
}

/*
 * IDALapackGetNumJacEvals returns the number of Jacobian evaluations.
 */
int IDALapackGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackGetNumJacEvals", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_LMEM_NULL, "IDALAPACK", "IDALapackGetNumJacEvals", MSGLS_LMEM_NULL);
    return(IDALAPACK_LMEM_NULL);
  }
  idalapack_mem = (IDALapackMem) lmem;

  *njevals = nje;

  return(IDALAPACK_SUCCESS);
}

/*
 * IDALapackGetNumResEvals returns the number of calls to the DAE function
 * needed for the DQ Jacobian approximation.
 */
int IDALapackGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackGetNumFctEvals", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_LMEM_NULL, "IDALAPACK", "IDALapackGetNumFctEvals", MSGLS_LMEM_NULL);
    return(IDALAPACK_LMEM_NULL);
  }
  idalapack_mem = (IDALapackMem) lmem;

  *nrevalsLS = nreDQ;

  return(IDALAPACK_SUCCESS);
}

/*
 * IDALapackGetReturnFlagName returns the name associated with a IDALAPACK
 * return value.
 */
char *IDALapackGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDALAPACK_SUCCESS:
    sprintf(name,"IDALAPACK_SUCCESS");
    break;   
  case IDALAPACK_MEM_NULL:
    sprintf(name,"IDALAPACK_MEM_NULL");
    break;
  case IDALAPACK_LMEM_NULL:
    sprintf(name,"IDALAPACK_LMEM_NULL");
    break;
  case IDALAPACK_ILL_INPUT:
    sprintf(name,"IDALAPACK_ILL_INPUT");
    break;
  case IDALAPACK_MEM_FAIL:
    sprintf(name,"IDALAPACK_MEM_FAIL");
    break;
  case IDALAPACK_JACFUNC_UNRECVR:
    sprintf(name,"IDALAPACK_JACFUNC_UNRECVR");
    break;
  case IDALAPACK_JACFUNC_RECVR:
    sprintf(name,"IDALAPACK_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * IDALapackGetLastFlag returns the last flag set in a IDALAPACK function.
 */
int IDALapackGetLastFlag(void *ida_mem, int *flag)
{
  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDALAPACK_MEM_NULL, "IDALAPACK", "IDALapackGetLastFlag", MSGLS_IDAMEM_NULL);
    return(IDALAPACK_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDALAPACK_LMEM_NULL, "IDALAPACK", "IDALapackGetLastFlag", MSGLS_LMEM_NULL);
    return(IDALAPACK_LMEM_NULL);
  }
  idalapack_mem = (IDALapackMem) lmem;

  *flag = last_flag;

  return(IDALAPACK_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * idaLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int idaLapackDenseInit(IDAMem IDA_mem)
{
  IDALapackMem idalapack_mem;

  idalapack_mem = (IDALapackMem) lmem;
  
  nje   = 0;
  nreDQ = 0;
  
  if (djac == NULL) {
    djac = idaLapackDenseDQJac;
    J_data = IDA_mem;
  }

  last_flag = IDALAPACK_SUCCESS;
  return(0);
}

/*
 * idaLapackDenseSetup does the setup operations for the dense linear solver. 
 * It calls the Jacobian function to obtain the Newton matrix M = F_y + c_j*F_y', 
 * updates counters, and calls the dense LU factorization routine.
 */
static int idaLapackDenseSetup(IDAMem IDA_mem,
                              N_Vector yP, N_Vector ypP, N_Vector fctP,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  IDALapackMem idalapack_mem;
  realtype fact;
  int ier, retval, one = 1;

  idalapack_mem = (IDALapackMem) lmem;

  /* Call Jacobian function */
  nje++;
  retval = djac(n, tn, cj, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDALAPACK_JACFUNC_UNRECVR, "IDALAPACK", "idaLapackDenseSetup", MSGLS_JACFUNC_FAILED);
    last_flag = IDALAPACK_JACFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    last_flag = IDALAPACK_JACFUNC_RECVR;
    return(1);
  }
  
  /* Do LU factorization of M */
  dgetrf_f77(&n, &n, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * idaLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int idaLapackDenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  IDALapackMem idalapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  idalapack_mem = (IDALapackMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &n, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1); 
  if (ier > 0) return(1);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) {
    fact = TWO/(ONE + cjratio);
    dscal_f77(&n, &fact, bd, &one); 
  }

  last_flag = IDALAPACK_SUCCESS;
  return(0);
}

/*
 * idaLapackDenseFree frees memory specific to the dense linear solver.
 */
static int idaLapackDenseFree(IDAMem IDA_mem)
{
  IDALapackMem  idalapack_mem;

  idalapack_mem = (IDALapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  free(idalapack_mem); 
  idalapack_mem = NULL;

  return(0);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * idaLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int idaLapackBandInit(IDAMem IDA_mem)
{
  IDALapackMem idalapack_mem;

  idalapack_mem = (IDALapackMem) lmem;

  nje   = 0;
  nreDQ = 0;

  if (bjac == NULL) {
    bjac = idaLapackBandDQJac;
    J_data = IDA_mem;
  }

  last_flag = IDALAPACK_SUCCESS;
  return(0);
}

/*
 * idaLapackBandSetup does the setup operations for the band linear solver.
 * It calls the Jacobian function to obtain the Newton matrix M = F_y + c_j*F_y', 
 * updates counters, and calls the band LU factorization routine.
 */
static int idaLapackBandSetup(IDAMem IDA_mem,
                             N_Vector yP, N_Vector ypP, N_Vector fctP, 
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  IDALapackMem idalapack_mem;
  realtype fact;
  int ier, retval, one = 1;

  idalapack_mem = (IDALapackMem) lmem;

  /* Call Jacobian function */
  nje++;
  retval = bjac(n, mu, ml, tn, cj, yP, ypP, fctP, M, J_data, tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDALAPACK_JACFUNC_UNRECVR, "IDALAPACK", "idaLapackBandSetup", MSGLS_JACFUNC_FAILED);
    last_flag = IDALAPACK_JACFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    last_flag = IDALAPACK_JACFUNC_RECVR;
    return(+1);
  }
  
  /* Do LU factorization of M */
  dgbtrf_f77(&n, &n, &ml, &mu, M->data, &(M->ldim), pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * idaLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int idaLapackBandSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  IDALapackMem idalapack_mem;
  realtype *bd, fact;
  int ier, one = 1;

  idalapack_mem = (IDALapackMem) lmem;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &n, &ml, &mu, &one, M->data, &(M->ldim), pivots, bd, &n, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in cj */
  if (cjratio != ONE) {
    fact = TWO/(ONE + cjratio);
    dscal_f77(&n, &fact, bd, &one); 
  }

  last_flag = IDALAPACK_SUCCESS;
  return(0);
}

/*
 * idaLapackBandFree frees memory specific to the band linear solver.
 */
static int idaLapackBandFree(IDAMem IDA_mem)
{
  IDALapackMem  idalapack_mem;

  idalapack_mem = (IDALapackMem) lmem;
  
  LapackFreeMat(M);
  LapackFreeArray(pivots);
  free(idalapack_mem); 
  idalapack_mem = NULL;

  return(0);
}

/* 
 * =================================================================
 *  DENSE DQ JACOBIAN APPROXIMATIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * idaLapackDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian F_y + c_j*F_y'. It assumes that a dense matrix of type
 * LapackMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro LAPACK_DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 
static int idaLapackDenseDQJac(int N, realtype tt, realtype c_j,
                               N_Vector yy, N_Vector yp, N_Vector rr, 
                               LapackMat Jac, void *jac_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  int j;
  int retval = 0;

  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* jac_data points to IDA_mem */
  IDA_mem = (IDAMem) jac_data;
  idalapack_mem = (IDALapackMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  rtemp  = tmp1;
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);
  if(constraints!=NULL) cns_data = N_VGetArrayPointer(constraints);

  srur = RSqrt(uround);

  for (j=0; j < N; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(LAPACK_DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as hh*yp_j. */

    inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewt_data[j] );

    if (hh*ypj < ZERO) inc = -inc;
    inc = (yj + inc) - yj;

    /* Adjust sign(inc) again if y_j has an inequality constraint. */
    if (constraints != NULL) {
      conj = cns_data[j];
      if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
      else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
    }

    /* Increment y_j and yp_j, call res, and break on error return. */
    y_data[j] += inc;
    yp_data[j] += c_j*inc;

    retval = res(tt, yy, yp, rtemp, rdata);
    nreDQ++;
    if (retval != 0) break;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, rr, jthCol);

    LAPACK_DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

    /*  reset y_j, yp_j */     
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

static int idaLapackBandDQJac(int N, int mupper, int mlower,
                              realtype tt, realtype c_j, 
                              N_Vector yy, N_Vector yp, N_Vector rr,
                              LapackMat Jac, void *jac_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj, ewtj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  realtype *ytemp_data, *yptemp_data, *rtemp_data, *r_data, *col_j;
  int group;
  
  N_Vector rtemp, ytemp, yptemp;
  int i, j, i1, i2, width, ngroups;
  int retval = 0;

  IDAMem IDA_mem;
  IDALapackMem idalapack_mem;

  /* jac_data points to IDA_mem */
  IDA_mem = (IDAMem) jac_data;
  idalapack_mem = (IDALapackMem) lmem;

  rtemp = tmp1; /* Rename work vector for use as the perturbed residual. */

  ytemp = tmp2; /* Rename work vector for use as a temporary for yy. */


  yptemp= tmp3; /* Rename work vector for use as a temporary for yp. */

  /* Obtain pointers to the data for all eight vectors used.  */

  ewt_data = N_VGetArrayPointer(ewt);
  r_data   = N_VGetArrayPointer(rr);
  y_data   = N_VGetArrayPointer(yy);
  yp_data  = N_VGetArrayPointer(yp);

  rtemp_data  = N_VGetArrayPointer(rtemp);
  ytemp_data  = N_VGetArrayPointer(ytemp);
  yptemp_data = N_VGetArrayPointer(yptemp);

  if (constraints != NULL) cns_data = N_VGetArrayPointer(constraints);

  /* Initialize ytemp and yptemp. */

  N_VScale(ONE, yy, ytemp);
  N_VScale(ONE, yp, yptemp);

  /* Compute miscellaneous values for the Jacobian computation. */

  srur = RSqrt(uround);
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all yy[j] and yp[j] for j in this group. */

    for (j=group-1; j<N; j+=width) {
        yj = y_data[j];
        ypj = yp_data[j];
        ewtj = ewt_data[j];

        /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
        adjustments using ypj and ewtj if this is small, and a further
        adjustment to give it the same sign as hh*ypj. */

        inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewtj );

        if (hh*ypj < ZERO) inc = -inc;
        inc = (yj + inc) - yj;

        /* Adjust sign(inc) again if yj has an inequality constraint. */

        if (constraints != NULL) {
          conj = cns_data[j];
          if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
          else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
        }

        /* Increment yj and ypj. */

        ytemp_data[j] += inc;
        yptemp_data[j] += cj*inc;
    }

    /* Call res routine with incremented arguments. */

    retval = res(tt, ytemp, yptemp, rtemp, rdata);
    nreDQ++;
    if (retval != 0) break;

    /* Loop over the indices j in this group again. */

    for (j=group-1; j<N; j+=width) {

      /* Reset ytemp and yptemp components that were perturbed. */

      yj = ytemp_data[j]  = y_data[j];
      ypj = yptemp_data[j] = yp_data[j];
      col_j = LAPACK_BAND_COL(Jac, j);
      ewtj = ewt_data[j];
      
      /* Set increment inc exactly as above. */

      inc = MAX( srur * MAX( ABS(yj), ABS(hh*ypj) ) , ONE/ewtj );
      if (hh*ypj < ZERO) inc = -inc;
      inc = (yj + inc) - yj;
      if (constraints != NULL) {
        conj = cns_data[j];
        if (ABS(conj) == ONE)      {if((yj+inc)*conj <  ZERO) inc = -inc;}
        else if (ABS(conj) == TWO) {if((yj+inc)*conj <= ZERO) inc = -inc;}
      }
      
      /* Load the difference quotient Jacobian elements for column j. */

      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower,N-1);
      
      for (i=i1; i<=i2; i++) 
            LAPACK_BAND_COL_ELEM(col_j,i,j) = inc_inv*(rtemp_data[i]-r_data[i]);
    }
    
  }
  
  return(retval);
  
}
