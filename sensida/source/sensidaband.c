/******************************************************************
 *                                                                *
 * File          : sensidaband.c                                  *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 3 July 2001                                    *
 *----------------------------------------------------------------*
 * This is the implementation file for the IDA banded linear      *
 * solver module, SENSIDABAND.                                    *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "sensidaband.h"
#include "ida.h"
#include "band.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "sensida.h"


/* Error Messages */

#define MSG_BAD_SIZES_1  "SensIDABand-- Illegal bandwidth parameter(s) "
#define MSG_BAD_SIZES_2  "mlower = %ld, mupper = %ld.\n"
#define MSG_BAD_SIZES_3  "Must have 0 <=  mlower, mupper <= N-1 =%ld.\n\n"
#define MSG_BAD_SIZES    MSG_BAD_SIZES_1 MSG_BAD_SIZES_2 MSG_BAD_SIZES_3

#define MSG_MEM_FAIL     "SensIDABandInit-- A memory request failed.\n\n"


/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/******************************************************************
 *                                                                *           
 * Types : IDABandMemRec, IDABandMem                              *
 *----------------------------------------------------------------*
 * The type IDABandMem is pointer to an IDABandMemRec. This       *
 * structure contains IDABand solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {

   IDABandJacFn b_jac;   /* jac = banded Jacobian routine to be called   */

   BandMat b_J;          /* J = dF/dy + cj*dF/dy', banded approximation. */
  
   integer b_mupper;     /* mupper = upper bandwidth of Jacobian matrix. */

   integer b_mlower;     /* mlower = lower bandwidth of Jacobian matrix. */

   integer b_storage_mu; /* storage_mu = upper bandwidth with storage for
                            factoring = min(N-1, mupper+mlower).       */

   integer *b_pivots;    /* pivots = pivot array for PJ = LU        */

   long int b_nje;       /* nje = no. of calls to jac               */
 
   void *b_jdata;        /* jdata = data structure required by jac. */

} IDABandMemRec, *IDABandMem;


/* SensIDABAND linit, lsetup, lsolve, and lfree routines */
 
static int  SensIDABandInit(IDAMem ida_mem, boole *setupNonNull);

static int  SensIDABandSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
			     N_Vector resp, N_Vector tempv1,
			     N_Vector tempv2, N_Vector tempv3);

static int  SensIDABandSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
			     N_Vector ypcur, N_Vector rescur);

static int SensIDABandFree(IDAMem ida_mem);

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

#define jac         (idaband_mem->b_jac)
#define JJ          (idaband_mem->b_J)
#define storage_mu  (idaband_mem->b_storage_mu)
#define pivots      (idaband_mem->b_pivots)
#define nje         (idaband_mem->b_nje)
#define jacdata     (idaband_mem->b_jdata)
 
                  
/*************** SensIDABand *****************************************

 This routine initializes the memory record and sets various function
 fields specific to the SENSIDABAND linear solver module.  SensIDABand sets
 the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and ida_lfree fields
 in (*IDA_mem) to be SensIDABandInit, SensIDABandSetup, SensIDABandSolve, 
 NULL, and SensIDABandFree, respectively.
 It allocates memory for a structure of type IDABandMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.
 It sets b_jdata field in the IDABandMemRec structure to be
 the input parameter jdata and the b_jac field to be:
   (1) the input parameter bjac, if bjac != NULL, or
   (2) IDABandDQJac, if bjac == NULL.
 Finally, it initializes SensIDABAND-specific optional outputs.

 IDABand returns SUCCESS = 0 if successful, or 
    SensIDA_BAND_FAIL    if IDA_mem = NULL or malloc failed, or
    SensIDA_BAND_BAD_ARG if mupper or mlower is illegal.
**********************************************************************/

int SensIDABand(void *IDA_mem, integer mupper, integer mlower,
               IDABandJacFn bjac, void *jdata)
{
  IDAMem ida_mem;
  IDABandMem idaband_mem;
  int N;
  SensData sdata;

  /* Return immediately if IDA_mem is NULL. */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) return(SensIDA_BAND_FAIL);

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  /* Set five main function fields in ida_mem. */
  linit  = SensIDABandInit;
  lsetup = SensIDABandSetup;
  lsolve = SensIDABandSolve;
  lperf  = NULL;
  lfree  = SensIDABandFree;

  /* Get memory for IDABandMemRec. */
  lmem = idaband_mem = (IDABandMem) malloc(sizeof(IDABandMemRec));
  if (idaband_mem == NULL) return(SensIDA_BAND_FAIL);

  /* Set Jacobian routine field to user's bjac or IDABandDQJac. */
  if (bjac == NULL) jac = IDABandDQJac;
  else jac = bjac;

  /* Test mlower and mupper for legality and load in memory. */
  if ((mlower < 0) || (mupper < 0) || (mlower >= N) || (mupper >= N)) {
    fprintf(errfp, MSG_BAD_SIZES, mlower, mupper, N-1);
    return(SensIDA_BAND_BAD_ARG);
  }
  idaband_mem->b_mlower = mlower;
  idaband_mem->b_mupper = mupper;
  jacdata = jdata;
     
  /* Set extended upper half-bandwidth for JJ (required for pivoting). */
  storage_mu = MIN(N-1, mupper + mlower);

  /* Initialize nje and set workspace lengths. */
  nje = 0;
  if (iopt != NULL) {
   iopt[BAND_NJE] = nje;
   iopt[BAND_LRW] = N*(storage_mu + mlower + 1);
   iopt[BAND_LIW] = N;
  }

  return(SUCCESS);
}

/* More convenience definitions */
#define mlower      (idaband_mem->b_mlower)
#define mupper      (idaband_mem->b_mupper)


/*************** SensIDABandInit *************************************

 This routine initializes remaining memory specific to the SENSIDABAND
 linear solver module.  If any memory request fails, memory previously
 allocated is freed, and an error message printed, before returning.
 The return value is either LINIT_OK (success) or LINIT_ERR (failure).

**********************************************************************/

static int SensIDABandInit(IDAMem ida_mem, boole *setupNonNull)
{
  IDABandMem idaband_mem;
  int N;
  SensData sdata;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  idaband_mem = (IDABandMem) lmem;

  /* Print error message and return if idaband_mem is NULL. */
  if (idaband_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }

  /* Set flag setupNonNull = TRUE. */
  *setupNonNull = TRUE;

  /* Allocate memory for JJ and pivot array. */
  
  JJ = BandAllocMat(N, mupper, mlower, storage_mu);
  if (JJ == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }
  pivots = BandAllocPiv(N);
  if (pivots == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    BandFreeMat(JJ);
    return(LINIT_ERR);
  }

   return(LINIT_OK);
}



/*************** SensIDABandSetup ****************************************

 This routine does the setup operations for the SENSIDABAND linear 
 solver module.  It calls the Jacobian evaluation routine,
 updates counters, and calls the band LU factorization routine.
 The return value is either
     SUCCESS = 0        if successful,
     LSETUP_ERROR_RECVR if the jac routine failed recoverably or the
                        LU factorization failed, or
     LSETUP_ERROR_NONRECVR if the jac routine failed unrecoverably.
**********************************************************************/

static int SensIDABandSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
			    N_Vector resp, N_Vector tempv1, N_Vector tempv2,
			    N_Vector tempv3)
{
  int retval;
  integer retfac;
  IDABandMem idaband_mem;
  int N;
  SensData sdata;
  N_Vector *yypsub, *yppsub, *respsub, *ewtsub;
  N_Vector *constraintssub, pconstraints;
  N_Vector *tempv1sub, *tempv2sub, *tempv3sub;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  idaband_mem = (IDABandMem) lmem;

  /* Increment nje counter. */
  nje++;

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
  
  if (iopt != NULL) iopt[BAND_NJE] = nje;

  /* Zero out JJ; call Jacobian routine jac; return if it failed. */
  BandZero(JJ);
  retval = jac(N, mupper, mlower, tn, yypsub[0], yppsub[0], cj, pconstraints, 
	       sdata->res_ptr, sdata->res_data, jacdata, respsub[0], ewtsub[0], 
	       hh, uround, JJ, &nre, tempv1sub[0], tempv2sub[0], tempv3sub[0]);
  if (retval < 0) return(LSETUP_ERROR_NONRECVR);
  if (retval > 0) return(LSETUP_ERROR_RECVR);

  /* Do LU factorization of JJ; return success or fail flag. */
  retfac = BandFactor(JJ, pivots);

  if (retfac != SUCCESS) return(LSETUP_ERROR_RECVR);
  return(SUCCESS);
}



/*************** SensIDABandSolve ****************************************

 This routine handles the solve operation for the SENSIDABAND linear
 solver module.  It calls the band backsolve routine, scales the
 solution vector according to cjratio, then returns SUCCESS = 0.

**********************************************************************/

static int SensIDABandSolve(IDAMem ida_mem, N_Vector b, N_Vector ycur,
			    N_Vector ypcur, N_Vector rescur)
{
  IDABandMem idaband_mem;
  int i, Ns;
  N_Vector *bsub;
  SensData sdata;
  
  sdata = (SensData) ida_mem->ida_rdata;
  Ns = sdata->Ns;
  bsub = N_VSUB(b);

  idaband_mem = (IDABandMem) lmem;

  for (i = 0; i <= Ns; i++) {
    BandBacksolve(JJ, pivots, bsub[i]);

    /* Scale the correction to account for change in cj. */
    if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), bsub[i], bsub[i]);
  }

  return(SUCCESS);
}



/*************** SensIDABandFree *****************************************

 This routine frees memory specific to the IDABAND linear solver.

**********************************************************************/

static int SensIDABandFree(IDAMem ida_mem)
{
  IDABandMem idaband_mem;

  idaband_mem = (IDABandMem) lmem;
  
  BandFreeMat(JJ);
  BandFreePiv(pivots);
  free(lmem);
}
