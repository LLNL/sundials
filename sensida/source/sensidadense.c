/******************************************************************
 *                                                                *
 * File          : sensidadense.c                                 *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 3 July 2001                                    *
 *----------------------------------------------------------------*
 * This is the implementation file for the SENSIDA dense linear   *
 * solver module, SENSIDADENSE.                                   *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "sensidadense.h"
#include "ida.h"
#include "dense.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "sensida.h"


/* Error Messages */

#define SensIDADENSE_INIT  "SensIDADenseInit-- "
#define MSG_MEM_FAIL  SensIDADENSE_INIT "A memory request failed.\n\n"


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
 
static int  SensIDADenseInit(IDAMem ida_mem, boole *setupNonNull);

static int  SensIDADenseSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
			      N_Vector resp, N_Vector tempv1,
			      N_Vector tempv2, N_Vector tempv3);

static int  SensIDADenseSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
			      N_Vector ypcur, N_Vector rescur);

static int  SensIDADenseFree(IDAMem ida_mem);

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

#define jac         (idadense_mem->d_jac)
#define JJ          (idadense_mem->d_J)
#define pivots      (idadense_mem->d_pivots)
#define nje         (idadense_mem->d_nje)
#define jacdata     (idadense_mem->d_jdata)
 
                  
/*************** SensIDADense *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the SENSIDADENSE linear solver module.  SENSIDADense sets
 the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and ida_lfree fields
 in (*IDA_mem) to be SensIDADenseInit, SensIDADenseSetup, SensIDADenseSolve, 
 NULL, and SensIDADenseFree, respectively.
 It allocates memory for a structure of type IDADenseMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.
 It sets d_jdata field in the IDADenseMemRec structure to be
 the input parameter jdata and the d_jac field to be:
   (1) the input parameter djac, if djac != NULL, or                
   (2) IDADenseDQJac, if djac == NULL.                             
 Finally, it initializes the IDADENSE-specific counters.

 SensIDADense returns 0 if successful, or SensIDA_DENSE_FAIL if either 
 IDA_mem is NULL or the call to malloc failed.
**********************************************************************/

int SensIDADense(void *IDA_mem, IDADenseJacFn djac, void *jdata)
{
  IDAMem ida_mem;
  IDADenseMem idadense_mem;
  int N;
  SensData sdata;

  /* Return immediately if IDA_mem is NULL. */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) return(SensIDA_DENSE_FAIL);

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  /* Set five main function fields in ida_mem. */
  linit  = SensIDADenseInit;
  lsetup = SensIDADenseSetup;
  lsolve = SensIDADenseSolve;
  lperf  = NULL;
  lfree  = SensIDADenseFree;

  /* Get memory for IDADenseMemRec. */
  lmem = idadense_mem = (IDADenseMem) malloc(sizeof(IDADenseMemRec));
  if (idadense_mem == NULL) return(SensIDA_DENSE_FAIL);

  /* Set Jacobian routine field to user's djac or IDADenseDQJac. */
  if (djac == NULL) jac = IDADenseDQJac;
  else jac = djac;

  jacdata = jdata;

  /* Initialize nje and set workspace lengths. */
   
  nje = 0;
  if (iopt != NULL) {
   iopt[DENSE_NJE] = nje;
   iopt[DENSE_LRW] = N*N;
   iopt[DENSE_LIW] = N;
  }

  return(SUCCESS);
}

/*************** SensIDADenseInit *************************************

 This routine initializes remaining memory specific to the SENSIDADENSE
 linear solver module.  If any memory request fails, memory previously
 allocated is freed, and an error message printed, before returning.
 The return value is either LINIT_OK (success) or LINIT_ERR (failure).

**********************************************************************/

static int SensIDADenseInit(IDAMem ida_mem, boole *setupNonNull)
{
  IDADenseMem idadense_mem;
  int N;
  SensData sdata;
  
  idadense_mem = (IDADenseMem) lmem;

  /* Print error message and return if idadense_mem is NULL. */
  if (idadense_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  /* Set flag setupNonNull = TRUE. */
  *setupNonNull = TRUE;
  
  /* Allocate memory for JJ and pivot array. */
  
  JJ = DenseAllocMat(N);
  if (JJ == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }
  pivots = DenseAllocPiv(N);
  if (pivots == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    DenseFreeMat(JJ);
    return(LINIT_ERR);
  }
  
  return(LINIT_OK);
}

/*************** SensIDADenseSetup ************************************

 This routine does the setup operations for the SENSIDADENSE linear 
 solver module.  It calls the Jacobian evaluation routine,
 updates counters, and calls the dense LU factorization routine.
 The return value is either
     SUCCESS = 0        if successful,
     LSETUP_ERROR_RECVR if the jac routine failed recoverably or the
                        LU factorization failed, or
     LSETUP_ERROR_NONRECVR if the jac routine failed unrecoverably.
**********************************************************************/

static int SensIDADenseSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
			     N_Vector resp, N_Vector tempv1, N_Vector tempv2,
			     N_Vector tempv3)
{
  int retval;
  integer retfac;
  IDADenseMem idadense_mem;
  int N;
  SensData sdata;
  N_Vector *yypsub, *yppsub, *respsub, *ewtsub;
  N_Vector *constraintssub, pconstraints;
  N_Vector *tempv1sub, *tempv2sub, *tempv3sub;
  
  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  idadense_mem = (IDADenseMem) lmem;

  /* Increment nje counter. */
  nje++;
  if (iopt != NULL) iopt[DENSE_NJE] = nje;

  yypsub = N_VSUB(yyp);
  yppsub = N_VSUB(ypp);
  respsub = N_VSUB(resp);
  ewtsub = N_VSUB(ewt);
  if (constraints != NULL) {
    constraintssub = N_VSUB(constraints);
    pconstraints = constraintssub[0];
  } 
  else {
    pconstraints = NULL;
  }
  tempv1sub = N_VSUB(tempv1);
  tempv2sub = N_VSUB(tempv2);
  tempv3sub = N_VSUB(tempv3);

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  DenseZero(JJ);
  retval = jac(N, tn, yypsub[0], yppsub[0], cj, pconstraints, sdata->res_ptr, 
	       sdata->res_data, jacdata, respsub[0], ewtsub[0], hh, uround, 
	       JJ, &nre, tempv1sub[0], tempv2sub[0], tempv3sub[0]);
  if (retval < 0) return(LSETUP_ERROR_NONRECVR);
  if (retval > 0) return(LSETUP_ERROR_RECVR);

  /* Do LU factorization of JJ; return success or fail flag. */
  retfac = DenseFactor(JJ, pivots);

  if (retfac != SUCCESS) return(LSETUP_ERROR_RECVR);
  return(SUCCESS);
}

/*************** SensIDADenseSolve ************************************

 This routine handles the solve operation for the SENSIDADENSE linear
 solver module.  It calls the dense backsolve routine, scales the
 solution vector according to cjratio, then returns SUCCESS = 0.

**********************************************************************/

static int SensIDADenseSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
			     N_Vector ypcur, N_Vector rescur)
{
  IDADenseMem idadense_mem;
  int i, Ns;
  N_Vector *bsub;
  SensData sdata;

  sdata = (SensData) ida_mem->ida_rdata;
  Ns = sdata->Ns;
  bsub = N_VSUB(b);

  idadense_mem = (IDADenseMem) lmem;

  for (i = 0; i <= Ns; i++) {
    DenseBacksolve(JJ, pivots, bsub[i]);

    /* Scale the correction to account for change in cj. */
    if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), bsub[i], bsub[i]);
  }

  return(SUCCESS);
}

/*************** SensIDADenseFree *****************************************

 This routine frees memory specific to the SENSIDADENSE linear solver.

**********************************************************************/

static int SensIDADenseFree(IDAMem ida_mem)
{
  IDADenseMem idadense_mem;

  idadense_mem = (IDADenseMem) lmem;
  
  DenseFreeMat(JJ);
  DenseFreePiv(pivots);
  free(lmem);
}
