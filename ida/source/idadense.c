/*******************************************************************
 *                                                                 *
 * File          : idadense.c                                      *
 * Programmers   : Alan C. Hindmarsh, Allan G. Taylor, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 11 July 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/ida/LICENSE                           *
 *-----------------------------------------------------------------*
 * This is the implementation file for the IDA dense linear        *
 * solver module, IDADENSE.                                        *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "idadense.h"
#include "ida.h"
#include "dense.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"

/* Error Messages */

#define IDADENSE               "IDADense-- "

#define MSG_IDAMEM_NULL        IDADENSE "IDA memory is NULL.\n\n"

#define MSG_MEM_FAIL           IDADENSE "A memory request failed.\n\n"

#define MSG_WRONG_NVEC         IDADENSE "Incompatible NVECTOR implementation.\n\n"

#define MSG_SETGET_IDAMEM_NULL "IDADenseSet*/IDADenseGet*-- IDA memory is NULL. \n\n"

#define MSG_SETGET_LMEM_NULL   "IDADenseSet*/IDADenseGet*-- IDADENSE memory is NULL. \n\n"

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* IDADENSE linit, lsetup, lsolve, and lfree routines */
 
static int IDADenseInit(IDAMem IDA_mem);

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1,
                         N_Vector tempv2, N_Vector tempv3);

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector ycur,
                         N_Vector ypcur, N_Vector rescur);

static int IDADenseFree(IDAMem IDA_mem);

static int IDADenseDQJac(integertype Neq, realtype tt, N_Vector yy, N_Vector yp,
                         realtype c_j, void *jdata, N_Vector resvec, DenseMat Jac,
                         N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

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
#define errfp        (IDA_mem->ida_errfp)
#define iopt         (IDA_mem->ida_iopt)
#define linit        (IDA_mem->ida_linit)
#define lsetup       (IDA_mem->ida_lsetup)
#define lsolve       (IDA_mem->ida_lsolve)
#define lperf        (IDA_mem->ida_lperf)
#define lfree        (IDA_mem->ida_lfree)
#define lmem         (IDA_mem->ida_lmem)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define nvspec       (IDA_mem->ida_nvspec)

#define neq          (idadense_mem->d_neq)
#define jac          (idadense_mem->d_jac)
#define JJ           (idadense_mem->d_J)
#define pivots       (idadense_mem->d_pivots)
#define nje          (idadense_mem->d_nje)
#define nreD         (idadense_mem->d_nreD)
#define jacdata      (idadense_mem->d_jdata)
 
                  
/*************** IDADense *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the IDADENSE linear solver module.  
 IDADense first calls the existing lfree routine if this is not NULL.
 Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 ida_lfree fields in (*IDA_mem) to be IDADenseInit, IDADenseSetup,
 IDADenseSolve, NULL, and IDADenseFree, respectively.
 It allocates memory for a structure of type IDADenseMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.
 It sets setupNonNull in (*IDA_mem) to TRUE, sets the d_jdata field
 in the IDADenseMemRec structure to be the input parameter jdata,
 and sets the d_jac field to be:
   (1) the input parameter djac, if djac != NULL, or                
   (2) IDADenseDQJac, if djac == NULL.                             
 Finally, it allocates memory for JJ and pivots.
 The return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

 NOTE: The dense linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, IDADense will first 
       test for compatible a compatible N_Vector internal
       representation by checking (1) the vector specification
       ID tag and (2) that the functions N_VMake, N_VDispose,
       N_VGetData, and N_VSetData are implemented.

**********************************************************************/

int IDADense(void *ida_mem, integertype Neq)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;
  int flag;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_IDAMEM_NULL);
    return(LMEM_FAIL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if ((strcmp(nvspec->tag,"serial")) || 
      nvspec->ops->nvmake    == NULL || 
      nvspec->ops->nvdispose == NULL ||
      nvspec->ops->nvgetdata == NULL || 
      nvspec->ops->nvsetdata == NULL) {
    fprintf(errfp, MSG_WRONG_NVEC);
    return(LIN_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  /* Get memory for IDADenseMemRec. */
  idadense_mem = (IDADenseMem) malloc(sizeof(IDADenseMemRec));
  if (idadense_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }

  /* Set default Jacobian routine and Jacobian data */
  jac = IDADenseDQJac;
  jacdata = IDA_mem;

  setupNonNull = TRUE;

  /* Store problem size */
  neq = Neq;

  /* Allocate memory for JJ and pivot array. */
  JJ = DenseAllocMat(Neq);
  if (JJ == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }
  pivots = DenseAllocPiv(Neq);
  if (pivots == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    DenseFreeMat(JJ);
    return(LMEM_FAIL);
  }

  /* Attach linear solver memory to IDA memory */
  lmem = idadense_mem;

  return(SUCCESS);
}

/*************  IDADenseSetJacFn **************************************/

int IDADenseSetJacFn(void *ida_mem, IDADenseJacFn djac)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  jac = djac;

  return(SUCCESS);
}

/*************** IDADenseSetJacData ***********************************/

int IDADenseSetJacData(void *ida_mem, void *jac_data)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  jacdata = jac_data;

  return(SUCCESS);
}

/*************** IDADenseGetIntWorkSpace ******************************/

int IDADenseGetIntWorkSpace(void *ida_mem, long int *leniwD)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  *leniwD = neq;

  return(OKAY);
}

/*************** IDADenseGetRealWorkSpace *****************************/

int IDADenseGetRealWorkSpace(void *ida_mem, long int *lenrwD)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  *lenrwD = neq*neq;

  return(OKAY);
}

/*************** IDADenseGetNumJacEvals *******************************/

int IDADenseGetNumJacEvals(void *ida_mem, int *njevalsD)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  *njevalsD = nje;

  return(OKAY);
}

/*************** IDADenseGetNumResEvals *******************************/

int IDADenseGetNumResEvals(void *ida_mem, int *nrevalsD)
{
  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stdout, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idadense_mem = (IDADenseMem) lmem;

  *nrevalsD = nreD;

  return(OKAY);
}


/*************** IDADenseInit *****************************************

 This routine does remaining initializations specific to the IDADENSE
 linear solver module.  It returns LINIT_OK = 0.

**********************************************************************/

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
  
  return(LINIT_OK);
}

/*************** IDADenseSetup ****************************************

 This routine does the setup operations for the IDADENSE linear 
 solver module.  It calls the Jacobian evaluation routine,
 updates counters, and calls the dense LU factorization routine.
 The return value is either
     SUCCESS = 0        if successful,
     LSETUP_ERROR_RECVR if the jac routine failed recoverably or the
                        LU factorization failed, or
     LSETUP_ERROR_NONRECVR if the jac routine failed unrecoverably.
**********************************************************************/

static int IDADenseSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1, N_Vector tempv2,
                         N_Vector tempv3)
{
  int retval;
  integertype retfac;
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Increment nje counter. */
  nje++;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  DenseZero(JJ);
  retval = jac(neq, tn, yyp, ypp, cj, jacdata, resp, JJ, 
               tempv1, tempv2, tempv3);
  if (retval < 0) return(LSETUP_ERROR_NONRECVR);
  if (retval > 0) return(LSETUP_ERROR_RECVR);

  /* Do LU factorization of JJ; return success or fail flag. */
  retfac = DenseFactor(JJ, pivots);

  if (retfac != SUCCESS) return(LSETUP_ERROR_RECVR);
  return(SUCCESS);
}

/*************** IDADenseSolve ****************************************

 This routine handles the solve operation for the IDADENSE linear
 solver module.  It calls the dense backsolve routine, scales the
 solution vector according to cjratio, then returns SUCCESS = 0.

**********************************************************************/

static int IDADenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector ycur,
                         N_Vector ypcur, N_Vector rescur)
{
  IDADenseMem idadense_mem;
  realtype *bd;
  
  idadense_mem = (IDADenseMem) lmem;
  
  bd = N_VGetData(b);
  DenseBacksolve(JJ, pivots, bd);
  N_VSetData(bd, b);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  return(SUCCESS);
}

/*************** IDADenseFree *****************************************

 This routine frees memory specific to the IDADENSE linear solver.

**********************************************************************/

static int IDADenseFree(IDAMem IDA_mem)
{
  IDADenseMem idadense_mem;

  idadense_mem = (IDADenseMem) lmem;
  
  DenseFreeMat(JJ);
  DenseFreePiv(pivots);
  free(lmem);

  return(0);

}

/*************** IDADenseDQJac ****************************************

 This routine generates a dense difference quotient approximation JJ to
 the DAE system Jacobian J.  It assumes that a dense matrix of type
 DenseMat is stored column-wise, and that elements within each column
 are contiguous.  The address of the jth column of JJ is obtained via
 the macro DENSE_COL and an N_Vector jthCol with the jth column as the
 component array is created using N_VMake and N_VGetData.  The jth column
 of the Jacobian is constructed using a call to the res routine, and
 a call to N_VLinearSum.
 The return value is either SUCCESS = 0, or the nonzero value returned
 by the res routine, if any.

**********************************************************************/

static int IDADenseDQJac(integertype Neq, realtype tt, N_Vector yy, N_Vector yp,
                         realtype c_j, void *jdata, N_Vector resvec, DenseMat Jac,
                         N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
 
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  integertype j;
  int retval=0;

  IDAMem IDA_mem;
  IDADenseMem idadense_mem;

  /* jdata points to IDA_mem */
  IDA_mem = (IDAMem) jdata;
  idadense_mem = (IDADenseMem) lmem;

  rtemp = tempv1; /* Rename work vector for use as residual value. */

  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VGetData(ewt);
  y_data   = N_VGetData(yy);
  yp_data  = N_VGetData(yp);
  if(constraints!=NULL) cns_data = N_VGetData(constraints);

  srur = RSqrt(uround);

  jthCol = N_VMake(NULL, nvspec); /* j loop overwrites this data address */

  for (j=0; j < Neq; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetData(DENSE_COL(Jac,j), jthCol);
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

    N_VSetData(y_data, yy);
    N_VSetData(yp_data, yp);

    retval = res(tt, yy, yp, rtemp, rdata);
    nreD++;
    if (retval != SUCCESS) break;

    /* Construct difference quotient in jthCol */
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, resvec, jthCol);
    DENSE_COL(Jac,j) = N_VGetData(jthCol);

    /*  reset y_j, yp_j */     
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  N_VDispose(jthCol);
  return(retval);

}

