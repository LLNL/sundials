/******************************************************************
 *                                                                *
 * File          : idadense.c                                     *
 * Programmers   : Alan C. Hindmarsh and Allan G. Taylor          *
 * Version of    : 30 August 1999                                 *
 *----------------------------------------------------------------*
 * This is the implementation file for the IDA dense linear       *
 * solver module, IDADENSE.                                       *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "idadense.h"
#include "ida.h"
#include "dense.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"


/* Error Messages */

#define IDADENSE_INIT  "IDADenseInit-- "
#define MSG_MEM_FAIL  IDADENSE_INIT "A memory request failed.\n\n"


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

    IDADenseJacFn d_jac; /* jac = Jacobian routine to be called  */

    DenseMat d_J;        /* J = dF/dy + cj*dF/dy'                */

    integer *d_pivots;   /* pivots = pivot array for PJ = LU     */

    long int d_nje;      /* nje = no. of calls to jac            */

    void *d_jdata;       /* jdata is passed to jac               */

} IDADenseMemRec, *IDADenseMem;


/* IDADENSE linit, lsetup, lsolve, and lfree routines */
 
static int  IDADenseInit(IDAMem ida_mem, boole *setupNonNull);

static int  IDADenseSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
                          N_Vector resp, N_Vector tempv1,
                          N_Vector tempv2, N_Vector tempv3);

static int  IDADenseSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
                          N_Vector ypcur, N_Vector rescur);

static int IDADenseFree(IDAMem ida_mem);


/*************** IDADenseDQJac ****************************************

 This routine generates a dense difference quotient approximation JJ to
 the DAE system Jacobian J.  It assumes that a dense matrix of type
 DenseMat is stored column-wise, and that elements within each column
 are contiguous.  The address of the jth column of JJ is obtained via
 the macro DENSE_COL and an N_Vector jthCol with the jth column as the
 component array is created using N_VMAKE and N_VDATA.  The jth column
 of the Jacobian is constructed using a call to the res routine, and
 a call to N_VLinearSum.
 The return value is either SUCCESS = 0, or the nonzero value returned
 by the res routine, if any.
**********************************************************************/

int IDADenseDQJac(integer Neq, real tt, N_Vector yy, N_Vector yp, real cj,
                  N_Vector constraints, ResFn res, void *rdata, void *jdata,
                  N_Vector resvec, N_Vector ewt, real hh, real uround,
                  DenseMat JJ, long int *nrePtr, N_Vector tempv1,
                  N_Vector tempv2, N_Vector tempv3)
 
{
  real inc, inc_inv, yj, ypj, srur, conj;
  real *y_data, *yp_data, *ewt_data, *cns_data;
  N_Vector rtemp, jthCol;
  integer j;
  int retval;

  rtemp = tempv1; /* Rename work vector for use as residual value. */

  /* Obtain pointers to the data for ewt, yy, yp. */
  ewt_data = N_VDATA(ewt);
  y_data   = N_VDATA(yy);
  yp_data  = N_VDATA(yp);
  if(constraints!=NULL)cns_data = N_VDATA(constraints);

  srur = RSqrt(uround);
  N_VMAKE(jthCol, NULL, Neq);

  for (j=0; j < Neq; j++) {

    /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

    /* Set data address of jthCol, and save y_j and yp_j values. */
    N_VDATA(jthCol) = DENSE_COL(JJ,j);
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
    (*nrePtr)++;
    retval = res(Neq, tt, yy, yp, rtemp, rdata);
    if (retval != SUCCESS) break;

    /* Construct difference quotient in jthCol, and reset y_j, yp_j. */ 
    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, rtemp, -inc_inv, resvec, jthCol);
    y_data[j] = yj;
    yp_data[j] = ypj;
  }

  N_VDISPOSE(jthCol);
  return(retval);

}


/* Readability Replacements */

#define Neq         (ida_mem->ida_Neq)
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

#define jac         (idadense_mem->d_jac)
#define JJ          (idadense_mem->d_J)
#define pivots      (idadense_mem->d_pivots)
#define nje         (idadense_mem->d_nje)
#define jacdata     (idadense_mem->d_jdata)
 
                  
/*************** IDADense *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the IDADENSE linear solver module.  IDADense sets
 the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and ida_lfree fields
 in (*IDA_mem) to be IDADenseInit, IDADenseSetup, IDADenseSolve, NULL,
 and IDADenseFree, respectively.
 It allocates memory for a structure of type IDADenseMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.
 It sets d_jdata field in the IDADenseMemRec structure to be
 the input parameter jdata and the d_jac field to be:
   (1) the input parameter djac, if djac != NULL, or                
   (2) IDADenseDQJac, if djac == NULL.                             
 Finally, it initializes the IDADENSE-specific counters.

 IDADense returns 0 if successful, or IDA_DENSE_FAIL if either 
 IDA_mem is NULL or the call to malloc failed.
**********************************************************************/

int IDADense(void *IDA_mem, IDADenseJacFn djac, void *jdata)
{
  IDAMem ida_mem;
  IDADenseMem idadense_mem;

  /* Return immediately if IDA_mem is NULL. */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) return(IDA_DENSE_FAIL);

  /* Set five main function fields in ida_mem. */
  linit  = IDADenseInit;
  lsetup = IDADenseSetup;
  lsolve = IDADenseSolve;
  lperf  = NULL;
  lfree  = IDADenseFree;

  /* Get memory for IDADenseMemRec. */
  lmem = idadense_mem = (IDADenseMem) malloc(sizeof(IDADenseMemRec));
  if (idadense_mem == NULL) return(IDA_DENSE_FAIL);

  /* Set Jacobian routine field to user's djac or IDADenseDQJac. */
  if (djac == NULL) jac = IDADenseDQJac;
  else jac = djac;

  jacdata = jdata;

  /* Initialize nje and set workspace lengths. */
   
  nje = 0;
  if (iopt != NULL) {
   iopt[DENSE_NJE] = nje;
   iopt[DENSE_LRW] = Neq*Neq;
   iopt[DENSE_LIW] = Neq;
  }

  return(SUCCESS);
}

/*************** IDADenseInit *****************************************

 This routine initializes remaining memory specific to the IDADENSE
 linear solver module.  If any memory request fails, memory previously
 allocated is freed, and an error message printed, before returning.
 The return value is either LINIT_OK (success) or LINIT_ERR (failure).

**********************************************************************/

static int IDADenseInit(IDAMem ida_mem, boole *setupNonNull)
{
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Print error message and return if idadense_mem is NULL. */
  if (idadense_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }

  /* Set flag setupNonNull = TRUE. */
  *setupNonNull = TRUE;
  
  /* Allocate memory for JJ and pivot array. */
  
  JJ = DenseAllocMat(Neq);
  if (JJ == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }
  pivots = DenseAllocPiv(Neq);
  if (pivots == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    DenseFreeMat(JJ);
    return(LINIT_ERR);
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
  integer retfac;
  IDADenseMem idadense_mem;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Increment nje counter. */
  nje++;
  if (iopt != NULL) iopt[DENSE_NJE] = nje;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  DenseZero(JJ);
  retval = jac(Neq, tn, yyp, ypp, cj, constraints, res, rdata, jacdata,
               resp, ewt, hh, uround, JJ, &nre, tempv1, tempv2, tempv3);
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
  
  idadense_mem = (IDADenseMem) lmem;
  
  DenseBacksolve(JJ, pivots, b);

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
}
