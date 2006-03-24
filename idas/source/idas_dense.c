/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-03-24 15:57:25 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the IDAS dense linear
 * solver module, IDADENSE.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "idas_dense_impl.h"
#include "idas_impl.h"
#include "idaa_impl.h"

#include "sundials_math.h"

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* IDADENSE linit, lsetup, lsolve, and lfree routines */
 
static int IDADenseInit(IDAMem IDA_mem);

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector rrp, N_Vector tmp1,
                         N_Vector tmp2, N_Vector tmp3);

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                         N_Vector ycur, N_Vector ypcur, N_Vector rrcur);

static int IDADenseFree(IDAMem IDA_mem);

/* IDADENSE lfreeB function */

static void IDADenseFreeB(IDAadjMem IDAADJ_mem);

/* Wrapper function for adjoint code */

static int IDAAdenseJac(long int NeqB, 
                        realtype tt, 
                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                        realtype c_jB, void *idaadj_mem, 
                        DenseMat JacB, 
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/* CVDENSE DQJac routine */

static int IDADenseDQJac(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
                         N_Vector rr, realtype c_j, void *jac_data, DenseMat Jac,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Readability Replacements */

#define res          (IDA_mem->ida_res)
#define rdata        (IDA_mem->ida_rdata)
#define uround       (IDA_mem->ida_uround)
#define tn           (IDA_mem->ida_tn)
#define hh           (IDA_mem->ida_hh)
#define cj           (IDA_mem->ida_cj)
#define cjratio      (IDA_mem->ida_cjratio)
#define ewt          (IDA_mem->ida_ewt)
#define constraints  (IDA_mem->ida_constraints)
#define nre          (IDA_mem->ida_nre)
#define iopt         (IDA_mem->ida_iopt)
#define linit        (IDA_mem->ida_linit)
#define lsetup       (IDA_mem->ida_lsetup)
#define lsolve       (IDA_mem->ida_lsolve)
#define lperf        (IDA_mem->ida_lperf)
#define lfree        (IDA_mem->ida_lfree)
#define lmem         (IDA_mem->ida_lmem)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define vec_tmpl     (IDA_mem->ida_tempv1)

#define neq          (idadense_mem->d_neq)
#define jac          (idadense_mem->d_jac)
#define JJ           (idadense_mem->d_J)
#define pivots       (idadense_mem->d_pivots)
#define nje          (idadense_mem->d_nje)
#define nreD         (idadense_mem->d_nreD)
#define jacdata      (idadense_mem->d_jdata)
#define last_flag    (idadense_mem->d_last_flag)
                  
/*
 * -----------------------------------------------------------------
 * IDADense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDADENSE linear solver module.  
 * IDADense first calls the existing lfree routine if this is not NULL.
 * Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDADenseInit, IDADenseSetup,
 * IDADenseSolve, NULL, and IDADenseFree, respectively.
 * It allocates memory for a structure of type IDADenseMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem) to TRUE, sets the d_jdata field
 * in the IDADenseMemRec structure to be the input parameter jdata,
 * and sets the d_jac field to be:
 *   (1) the input parameter djac, if djac != NULL, or                
 *   (2) IDADenseDQJac, if djac == NULL.                             
 * Finally, it allocates memory for JJ and pivots.
 * The return value is IDADENSE_SUCCESS = 0, IDADENSE_LMEM_FAIL = -1,
 * or IDADENSE_ILL_INPUT = -2.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDADense will first 
 *       test for a compatible N_Vector internal
 *       representation by checking that the functions N_VGetArrayPointer
 *       and N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int IDADense(void *ida_mem, long int Neq)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;
  int flag;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADense", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if(vec_tmpl->ops->nvgetarraypointer == NULL ||
     vec_tmpl->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_ILL_INPUT, "IDADENSE", "IDADense", MSGD_BAD_NVECTOR);
    return(IDADENSE_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  /* Get memory for IDADenseMemRec. */
  idadense_mem = NULL;
  idadense_mem = (IDADenseMem) malloc(sizeof(IDADenseMemRec));
  if (idadense_mem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    return(IDADENSE_MEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  jac = IDADenseDQJac;
  jacdata = IDA_mem;
  last_flag = IDADENSE_SUCCESS;

  setupNonNull = TRUE;

  /* Store problem size */
  neq = Neq;

  /* Allocate memory for JJ and pivot array. */
  JJ = NULL;
  JJ = DenseAllocMat(Neq);
  if (JJ == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    free(idadense_mem); idadense_mem = NULL;
    return(IDADENSE_MEM_FAIL);
  }

  pivots = NULL;
  pivots = DenseAllocPiv(Neq);
  if (pivots == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_MEM_FAIL, "IDADENSE", "IDADense", MSGD_MEM_FAIL);
    DenseFreeMat(JJ);
    free(idadense_mem); idadense_mem = NULL;
    return(IDADENSE_MEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  lmem = idadense_mem;

  return(IDADENSE_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDADenseSet* and IDADenseGet*
 * -----------------------------------------------------------------
 */

int IDADenseSetJacFn(void *ida_mem, IDADenseJacFn djac, void *jac_data)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADenseSetJacFn", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_LMEM_NULL, "IDADENSE", "IDADenseSetJacFn", MSGD_LMEM_NULL);
    return(IDADENSE_LMEM_NULL);
  }
  idadense_mem = (IDADenseMem) lmem;

  jac = djac;
  if (djac != NULL) jacdata = jac_data;

  return(IDADENSE_SUCCESS);
}

int IDADenseGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADenseGetWorkSpace", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_LMEM_NULL, "IDADENSE", "IDADenseGetWorkSpace", MSGD_LMEM_NULL);
    return(IDADENSE_LMEM_NULL);
  }
  idadense_mem = (IDADenseMem) lmem;

  *lenrwLS = neq*neq;
  *leniwLS = neq;
 
  return(IDADENSE_SUCCESS);
}

int IDADenseGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADenseGetNumJacEvals", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_LMEM_NULL, "IDADENSE", "IDADenseGetNumJacEvals", MSGD_LMEM_NULL);
    return(IDADENSE_LMEM_NULL);
  }
  idadense_mem = (IDADenseMem) lmem;

  *njevals = nje;

  return(IDADENSE_SUCCESS);
}

int IDADenseGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADenseGetNumResEvals", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_LMEM_NULL, "IDADENSE", "IDADenseGetNumResEvals", MSGD_LMEM_NULL);
    return(IDADENSE_LMEM_NULL);
  }
  idadense_mem = (IDADenseMem) lmem;

  *nrevalsLS = nreD;

  return(IDADENSE_SUCCESS);
}

int IDADenseGetLastFlag(void *ida_mem, int *flag)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_MEM_NULL, "IDADENSE", "IDADenseGetLastFlag", MSGD_IDAMEM_NULL);
    return(IDADENSE_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    IDAProcessError(IDA_mem, IDADENSE_LMEM_NULL, "IDADENSE", "IDADenseGetLastFlag", MSGD_LMEM_NULL);
    return(IDADENSE_LMEM_NULL);
  }
  idadense_mem = (IDADenseMem) lmem;

  *flag = last_flag;

  return(IDADENSE_SUCCESS);
}

char *IDADenseGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDADENSE_SUCCESS:
    sprintf(name,"IDADENSE_SUCCESS");
    break;   
  case IDADENSE_MEM_NULL:
    sprintf(name,"IDADENSE_MEM_NULL");
    break;
  case IDADENSE_LMEM_NULL:
    sprintf(name,"IDADENSE_LMEM_NULL");
    break;
  case IDADENSE_ILL_INPUT:
    sprintf(name,"IDADENSE_ILL_INPUT");
    break;
  case IDADENSE_MEM_FAIL:
    sprintf(name,"IDADENSE_MEM_FAIL");
    break;
  case IDADENSE_JACFUNC_UNRECVR:
    sprintf(name,"IDADENSE_JACFUNC_UNRECVR");
    break;
  case IDADENSE_JACFUNC_RECVR:
    sprintf(name,"IDADENSE_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * IDADENSE interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the IDADENSE
  linear solver module.  It returns 0.
*/

static int IDADenseInit(IDAMem IDA_mem)
{
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

   
  nje  = 0;
  nreD = 0;

  if (jac == NULL) {
    jac = IDADenseDQJac;
    jacdata = IDA_mem;
  }
  
  last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the IDADENSE linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the dense LU factorization routine.
  The return value is either
     IDADENSE_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector rrp, N_Vector tmp1, N_Vector tmp2,
                         N_Vector tmp3)
{
  int retval;
  long int retfac;
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Increment nje counter. */
  nje++;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  DenseZero(JJ);
  retval = jac(neq, tn, yyp, ypp, rrp, cj, jacdata, JJ, 
               tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADENSE_JACFUNC_UNRECVR, "IDADENSE", "IDADenseSetup", MSGD_JACFUNC_FAILED);
    last_flag = IDADENSE_JACFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    last_flag = IDADENSE_JACFUNC_RECVR;
    return(+1);
  }

  /* Do LU factorization of JJ; return success or fail flag. */
  retfac = DenseFactor(JJ, pivots);

  if (retfac != 0) {
    last_flag = retfac;
    return(+1);
  }
  last_flag = IDADENSE_SUCCESS;
  return(0);
}

/*
  This routine handles the solve operation for the IDADENSE linear
  solver module.  It calls the dense backsolve routine, scales the
  solution vector according to cjratio, then returns IDADENSE_SUCCESS = 0.
*/

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                         N_Vector ycur, N_Vector ypcur, N_Vector rrcur)
{
  IDADenseMem idadense_mem;
  realtype *bd;
  
  idadense_mem = (IDADenseMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseBacksolve(JJ, pivots, bd);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  last_flag = 0;
  return(0);
}

/*
  This routine frees memory specific to the IDADENSE linear solver.
*/

static int IDADenseFree(IDAMem IDA_mem)
{
  IDADenseMem idadense_mem;

  idadense_mem = (IDADenseMem) lmem;
  
  DenseFreeMat(JJ);
  DenseFreePiv(pivots);
  free(lmem); lmem = NULL;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * IDADENSE private routines
 * -----------------------------------------------------------------
 */

/*
  This routine generates a dense difference quotient approximation Jac to
  the DAE system Jacobian J.  It assumes that a dense matrix of type
  DenseMat is stored column-wise, and that elements within each column
  are contiguous.  The address of the jth column of Jac is obtained via
  the macro DENSE_COL and this pointer is associated with an N_Vector
  using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
  The jth column of the Jacobian is constructed using a call to the res 
  routine, and a call to N_VLinearSum.
  The return value is either IDADENSE_SUCCESS = 0, or the nonzero value returned
  by the res routine, if any.
*/

static int IDADenseDQJac(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
                         N_Vector rr, realtype c_j, void *jac_data, DenseMat Jac,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *tmp2_data, *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  long int j;
  int retval=0;

  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* jac_data points to IDA_mem */
  IDA_mem = (IDAMem) jac_data;
  idadense_mem = (IDADenseMem) lmem;

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

  for (j=0; j < Neq; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);
    yj = y_data[j];
    ypj = yp_data[j];

    /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
    adjustments using yp_j and ewt_j if this is small, and a further
    adjustment to give it the same sign as hh*yp_j. */

    inc = srur*MAX(ABS(yj),MAX( ABS(hh*ypj), ONE/ewt_data[j]));
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
    nreD++;
    if (retval != IDADENSE_SUCCESS) break;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, rr, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);

    /*  reset y_j, yp_j */     
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);

}


/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/* Additional readability replacements */

#define ytmp        (IDAADJ_mem->ia_ytmp)
#define yptmp       (IDAADJ_mem->ia_yptmp)
#define getY        (IDAADJ_mem->ia_getY)
#define lmemB       (IDAADJ_mem->ia_lmemB)
#define lfreeB      (IDAADJ_mem->ia_lfreeB)

#define djac_B      (idadenseB_mem->d_djacB)
#define jac_data_B  (idadenseB_mem->d_jac_dataB)

/*
 * IDADenseB and IDADenseSet*B
 *
 * Wrappers for the backward phase around the corresponding 
 * IDAS functions
 */

int IDADenseB(void *idaadj_mem, long int NeqB)
{
  IDAadjMem IDAADJ_mem;
  IDADenseMemB idadenseB_mem;
  IDAMem IDAB_mem;
  int flag;

  if (idaadj_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_ADJMEM_NULL, "IDADENSE", "IDADenseB", MSGD_CAMEM_NULL);
    return(IDADENSE_ADJMEM_NULL);
  }
  IDAADJ_mem = (IDAadjMem) idaadj_mem;

  IDAB_mem = (IDAMem) IDAADJ_mem->IDAB_mem;

  /* Get memory for IDADenseMemRecB */
  idadenseB_mem = (IDADenseMemB) malloc(sizeof(IDADenseMemRecB));
  if (idadenseB_mem == NULL) {
    IDAProcessError(IDAB_mem, IDADENSE_MEM_FAIL, "IDADENSE", "IDADenseB", MSGD_MEM_FAIL);
    return(IDADENSE_MEM_FAIL);
  }

  djac_B = NULL;
  jac_data_B = NULL;

  /* attach lmemB and lfreeB */
  lmemB = idadenseB_mem;
  lfreeB = IDADenseFreeB;

  flag = IDADense(IDAB_mem, NeqB);

  if (flag != IDADENSE_SUCCESS) {
    free(idadenseB_mem);
    idadenseB_mem = NULL;
  }

  return(flag);
}

int IDADenseSetJacFnB(void *idaadj_mem, IDADenseJacFnB djacB, void *jdataB)
{
  IDAadjMem IDAADJ_mem;
  IDADenseMemB idadenseB_mem;
  IDAMem IDAB_mem;
  int flag;

  if (idaadj_mem == NULL) {
    IDAProcessError(NULL, IDADENSE_ADJMEM_NULL, "IDADENSE", "IDADenseSetJacFnB", MSGD_CAMEM_NULL);
    return(IDADENSE_ADJMEM_NULL);
  }
  IDAADJ_mem = (IDAadjMem) idaadj_mem;

  IDAB_mem = (IDAMem) IDAADJ_mem->IDAB_mem;

  if (lmemB == NULL) {
    IDAProcessError(IDAB_mem, IDADENSE_LMEMB_NULL, "IDADENSE", "IDADenseSetJacFnB", MSGD_LMEMB_NULL);
    return(IDADENSE_LMEMB_NULL);
  }
  idadenseB_mem = (IDADenseMemB) lmemB;

  djac_B     = djacB;
  jac_data_B = jdataB;

  flag = IDADenseSetJacFn(IDAB_mem, IDAAdenseJac, idaadj_mem);

  return(flag);
}

/*
 * IDADenseFreeB 
 */

static void IDADenseFreeB(IDAadjMem IDAADJ_mem)
{
  IDADenseMemB idadenseB_mem;

  idadenseB_mem = (IDADenseMemB) lmemB;

  free(idadenseB_mem);
}

/*
 * IDAAdenseJac
 *
 * This routine interfaces to the IDADenseJacFnB routine provided 
 * by the user.
 * NOTE: jac_data actually contains idaadj_mem
 */

static int IDAAdenseJac(long int NeqB, realtype tt, 
                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                        realtype c_jB, void *idaadj_mem, 
                        DenseMat JacB, 
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  IDAadjMem IDAADJ_mem;
  IDAMem IDAB_mem;
  IDADenseMemB idadenseB_mem;
  int flag;

  IDAADJ_mem = (IDAadjMem) idaadj_mem;
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  idadenseB_mem = (IDADenseMemB) lmemB;

  /* Forward solution from interpolation */
  flag = getY(IDAADJ_mem, tt, ytmp, yptmp);
  if (flag != IDA_SUCCESS) {
    IDAProcessError(IDAB_mem, -1, "IDADENSE", "IDAAdenseJac", MSGD_BAD_T);
    return(-1);
  }

  /* Call user's adjoint dense djacB routine */
  flag = djac_B(NeqB, tt, 
                ytmp, yptmp, 
                yyB, ypB, rrB, 
                c_jB, jac_data_B, 
                JacB, 
                tmp1B, tmp2B, tmp3B);

  return(flag);
}
