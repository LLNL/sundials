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

#define IDADENSE         "IDADense/IDAReInitDense-- "

#define MSG_IDAMEM_NULL  IDADENSE "IDA memory is NULL.\n\n"

#define MSG_MEM_FAIL     IDADENSE "A memory request failed.\n\n"

#define MSG_WRONG_NVEC   IDADENSE "Incompatible NVECTOR implementation.\n\n"

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)


/******************************************************************
 *                                                                *           
 * Types : IDADenseMemRec, IDADenseMem                            *
 *----------------------------------------------------------------*
 * The type IDADenseMem is pointer to an IDADenseMemRec. This     *
 * structure contains IDADense solver-specific data.              *
 *                                                                *
 ******************************************************************/

typedef struct {

  integertype d_Neq;     /* Neq = problem dimension              */

  IDADenseJacFn d_jac;   /* jac = Jacobian routine to be called  */
  
  DenseMat d_J;          /* J = dF/dy + cj*dF/dy'                */
  
  integertype *d_pivots; /* pivots = pivot array for PJ = LU     */
  
  long int d_nje;        /* nje = no. of calls to jac            */
  
  void *d_jdata;         /* jdata is passed to jac               */

} IDADenseMemRec, *IDADenseMem;


/* IDADENSE linit, lsetup, lsolve, and lfree routines */
 
static int IDADenseInit(IDAMem ida_mem);

static int IDADenseSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1,
                         N_Vector tempv2, N_Vector tempv3);

static int IDADenseSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
                         N_Vector ypcur, N_Vector rescur);

static int IDADenseFree(IDAMem ida_mem);

static int IDADenseDQJac(integertype neq, realtype tt, N_Vector yy, N_Vector yp,
                         realtype cj, void *jdata, N_Vector resvec, DenseMat JJ,
                         N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);



/* Readability Replacements */

#define res         (ida_mem->ida_res)
#define rdata       (ida_mem->ida_rdata)
#define uround      (ida_mem->ida_uround)
#define tn          (ida_mem->ida_tn)
#define hh          (ida_mem->ida_hh)
#define cj          (ida_mem->ida_cj)
#define cjratio     (ida_mem->ida_cjratio)
#define ewt         (ida_mem->ida_ewt)
#define constraints (ida_mem->ida_constraints)
#define nre         (ida_mem->ida_nre)
#define errfp       (ida_mem->ida_errfp)
#define iopt        (ida_mem->ida_iopt)
#define linit       (ida_mem->ida_linit)
#define lsetup      (ida_mem->ida_lsetup)
#define lsolve      (ida_mem->ida_lsolve)
#define lperf       (ida_mem->ida_lperf)
#define lfree       (ida_mem->ida_lfree)
#define lmem        (ida_mem->ida_lmem)
#define machenv     (ida_mem->ida_machenv)
#define setupNonNull  (ida_mem->ida_setupNonNull)

#define Neq         (idadense_mem->d_Neq)
#define jac         (idadense_mem->d_jac)
#define JJ          (idadense_mem->d_J)
#define pivots      (idadense_mem->d_pivots)
#define nje         (idadense_mem->d_nje)
#define jacdata     (idadense_mem->d_jdata)
 
                  
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
       representation by checking (1) the machine environment
       ID tag and (2) that the functions N_VMake, N_VDispose,
       N_VGetData, and N_VSetData are implemented.

**********************************************************************/

int IDADense(void *IDA_mem, integertype neq, IDADenseJacFn djac, void *jdata)
{
  IDAMem ida_mem;
  IDADenseMem idadense_mem;
  int flag;

  /* Return immediately if IDA_mem is NULL. */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) {
    fprintf(errfp, MSG_IDAMEM_NULL);
    return(LMEM_FAIL);
  }

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if ((strcmp(machenv->tag,"serial")) || 
      machenv->ops->nvmake    == NULL || 
      machenv->ops->nvdispose == NULL ||
      machenv->ops->nvgetdata == NULL || 
      machenv->ops->nvsetdata == NULL) {
    fprintf(errfp, MSG_WRONG_NVEC);
    return(LIN_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree(ida_mem);

  /* Set five main function fields in ida_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  /* Get memory for IDADenseMemRec. */
  lmem = idadense_mem = (IDADenseMem) malloc(sizeof(IDADenseMemRec));
  if (idadense_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }

  /* Set Jacobian routine field to user's djac or IDADenseDQJac. */
  if (djac == NULL) {
    jac = IDADenseDQJac;
    jacdata = IDA_mem;
  } else {
    jac = djac;
    jacdata = jdata;
  }
  setupNonNull = TRUE;

  /* Store problem size */
  Neq = neq;

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
  
  return(SUCCESS);
}

/*************** IDAReInitDense****************************************

 This routine resets the link between the main IDA module and the
 dense linear solver module IDADENSE.  No memory freeing or allocation
 operations are done, as the existing linear solver memory is assumed
 sufficient.  All other initializations are the same as in IDADense.
 The return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

**********************************************************************/

int IDAReInitDense(void *IDA_mem, IDADenseJacFn djac, void *jdata)
{
  IDAMem ida_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if IDA_mem is NULL. */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) {
    fprintf(errfp, MSG_IDAMEM_NULL);
    return(LMEM_FAIL);
  }

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if ((strcmp(machenv->tag,"serial")) || 
      machenv->ops->nvmake    == NULL || 
      machenv->ops->nvdispose == NULL ||
      machenv->ops->nvgetdata == NULL || 
      machenv->ops->nvsetdata == NULL) {
    fprintf(errfp, MSG_WRONG_NVEC);
    return(LIN_ILL_INPUT);
  }

  /* Set five main function fields in ida_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  idadense_mem = lmem;  /* Use existing linear solver memory pointer. */

  /* Set Jacobian routine field to user's djac or IDADenseDQJac. */
  if (djac == NULL) {
    jac = IDADenseDQJac;
    jacdata = IDA_mem;
  } else {
    jac = djac;
    jacdata = jdata;
  }
  setupNonNull = TRUE;

  return(SUCCESS);
}


/*************** IDADenseInit *****************************************

 This routine does remaining initializations specific to the IDADENSE
 linear solver module.  It returns LINIT_OK = 0.

**********************************************************************/

static int IDADenseInit(IDAMem ida_mem)
{
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Initialize nje and set workspace lengths. */
   
  nje = 0;
  if (iopt != NULL) {
   iopt[DENSE_NJE] = nje;
   iopt[DENSE_LRW] = Neq*Neq;
   iopt[DENSE_LIW] = Neq;
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

static int IDADenseSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1, N_Vector tempv2,
                         N_Vector tempv3)
{
  int retval;
  integertype retfac;
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Increment nje counter. */
  nje++;
  if (iopt != NULL) iopt[DENSE_NJE] = nje;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  DenseZero(JJ);
  retval = jac(Neq, tn, yyp, ypp, cj, jacdata, resp, JJ, 
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

static int IDADenseSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
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

static int IDADenseFree(IDAMem ida_mem)
{
  IDADenseMem idadense_mem;

  idadense_mem = (IDADenseMem) lmem;
  
  DenseFreeMat(JJ);
  DenseFreePiv(pivots);
  free(lmem);

  return(0);

}

#undef cj
#undef JJ


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

static int IDADenseDQJac(integertype neq, realtype tt, N_Vector yy, N_Vector yp,
                         realtype cj, void *jdata, N_Vector resvec, DenseMat JJ,
                         N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
 
{
  realtype inc, inc_inv, yj, ypj, srur, conj;
  realtype *y_data, *yp_data, *ewt_data, *cns_data = NULL;
  N_Vector rtemp, jthCol;
  integertype j;
  int retval;

  IDAMem ida_mem;

  /* jdata points to IDA_mem */
  ida_mem = (IDAMem) jdata;

  rtemp = tempv1; /* Rename work vector for use as residual value. */

  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VGetData(ewt);
  y_data   = N_VGetData(yy);
  yp_data  = N_VGetData(yp);
  if(constraints!=NULL) cns_data = N_VGetData(constraints);

  srur = RSqrt(uround);

  jthCol = N_VMake(NULL, machenv); /* j loop overwrites this data address */

  for (j=0; j < neq; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VSetData(DENSE_COL(JJ,j), jthCol);
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
    yp_data[j] += cj*inc;
    retval = res(tt, yy, yp, rtemp, rdata);
    nre++;
    if (retval != SUCCESS) break;

    /* Construct difference quotient in jthCol, and reset y_j, yp_j. */ 
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, resvec, jthCol);
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  N_VDispose(jthCol);
  return(retval);

}

