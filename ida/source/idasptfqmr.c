/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2005-10-26 23:08:08 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the IDA scaled preconditioned
 * TFQMR linear solver module, IDASPTFQMR.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "idasptfqmr_impl.h"
#include "sundialsmath.h"

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define PT9  RCONST(0.9)
#define PT05 RCONST(0.05)

#define IDA_SPTFQMR_MAXL 5

/* IDASPTFQMR linit, lsetup, lsolve, lperf, and lfree routines */

static int IDASptfqmrInit(IDAMem IDA_mem);

static int IDASptfqmrSetup(IDAMem IDA_mem, 
			   N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDASptfqmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
			   N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDASptfqmrPerf(IDAMem IDA_mem, int perftask);

static int IDASptfqmrFree(IDAMem IDA_mem);

/* IDASPTFQMR Atimes and PSolve routines called by generic SPTFQMR solver */

static int IDASptfqmrAtimes(void *ida_mem, N_Vector v, N_Vector z);

static int IDASptfqmrPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

static int IDASptfqmrDQJtimes(realtype tt,
			      N_Vector yy, N_Vector yp, N_Vector rr,
			      N_Vector v, N_Vector Jv, 
			      realtype c_j, void *jac_data, 
			      N_Vector work1, N_Vector work2);

/* Readability Replacements */

#define lrw1         (IDA_mem->ida_lrw1)
#define liw1         (IDA_mem->ida_liw1)
#define nst          (IDA_mem->ida_nst)
#define tn           (IDA_mem->ida_tn)
#define cj           (IDA_mem->ida_cj)
#define epsNewt      (IDA_mem->ida_epsNewt)
#define nre          (IDA_mem->ida_nre)
#define res          (IDA_mem->ida_res)
#define rdata        (IDA_mem->ida_rdata)
#define ewt          (IDA_mem->ida_ewt)
#define errfp        (IDA_mem->ida_errfp)
#define iopt         (IDA_mem->ida_iopt)
#define linit        (IDA_mem->ida_linit)
#define lsetup       (IDA_mem->ida_lsetup)
#define lsolve       (IDA_mem->ida_lsolve)
#define lperf        (IDA_mem->ida_lperf)
#define lfree        (IDA_mem->ida_lfree)
#define lmem         (IDA_mem->ida_lmem)
#define nni          (IDA_mem->ida_nni)
#define ncfn         (IDA_mem->ida_ncfn)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define vec_tmpl     (IDA_mem->ida_tempv1)

#define sqrtN       (idasptfqmr_mem->q_sqrtN)
#define epslin      (idasptfqmr_mem->q_epslin)
#define ytemp       (idasptfqmr_mem->q_ytemp)
#define yptemp      (idasptfqmr_mem->q_yptemp)
#define xx          (idasptfqmr_mem->q_xx)
#define ycur        (idasptfqmr_mem->q_ycur)
#define ypcur       (idasptfqmr_mem->q_ypcur)
#define rcur        (idasptfqmr_mem->q_rcur)
#define resflag     (idasptfqmr_mem->q_resflag)
#define npe         (idasptfqmr_mem->q_npe)
#define nli         (idasptfqmr_mem->q_nli)
#define nps         (idasptfqmr_mem->q_nps)
#define ncfl        (idasptfqmr_mem->q_ncfl)
#define nst0        (idasptfqmr_mem->q_nst0)
#define nni0        (idasptfqmr_mem->q_nni0)
#define nli0        (idasptfqmr_mem->q_nli0)
#define ncfn0       (idasptfqmr_mem->q_ncfn0)
#define ncfl0       (idasptfqmr_mem->q_ncfl0)
#define nwarn       (idasptfqmr_mem->q_nwarn)
#define njtimes     (idasptfqmr_mem->q_njtimes)
#define nreSQ       (idasptfqmr_mem->q_nreSQ)
#define sptfqmr_mem (idasptfqmr_mem->q_sptfqmr_mem)
#define last_flag   (idasptfqmr_mem->q_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmr
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDASPTFQMR linear solver module.
 *
 * IDASptfqmr first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDASptfqmrInit, IDASptfqmrSetup,
 * IDASptfqmrSolve, IDASptfqmrPerf, and IDASptfqmrFree, respectively.
 * It allocates memory for a structure of type IDASptfqmrMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem). It then sets the following
 * fields in the IDASptfqmrMemRec structure:
 *
 *   q_maxl      = IDA_SPTFQMR_MAXL  if maxl <= 0
 *               = maxl              if maxl > 0
 *   q_eplifac   = PT05
 *   q_dqincfac  = ONE
 *   q_pdata     = NULL
 *   q_pset      = NULL
 *   q_psolve    = NULL
 *   q_jtimes    = IDASptfqmrDQJtimes
 *   q_jdata     = ida_mem
 *   q_last_flag = IDASPTFQMR_SUCCESS
 *
 * Finally, IDASptfqmr allocates memory for ytemp, yptemp, and xx, and
 * calls SptfqmrMalloc to allocate memory for the Sptfqmr solver.
 *
 * The return value of IDASptfqmr is:
 *   IDASPTFQMR_SUCCESS   =  0 if successful
 *   IDASPTFQMR_MEM_FAIL  = -1 if IDA_mem is NULL or a memory
 *                             allocation failed
 *   IDASPTFQMR_ILL_INPUT = -2 if a required vector operation is not
 *                             implemented.
 * -----------------------------------------------------------------
 */

int IDASptfqmr(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;
  int flag, maxl1;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_BAD_NVECTOR);
    return(IDASPTFQMR_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  linit  = IDASptfqmrInit;
  lsetup = IDASptfqmrSetup;
  lsolve = IDASptfqmrSolve;
  lperf  = IDASptfqmrPerf;
  lfree  = IDASptfqmrFree;

  /* Get memory for IDASptfqmrMemRec */
  idasptfqmr_mem = (IDASptfqmrMem) malloc(sizeof(IDASptfqmrMemRec));
  if (idasptfqmr_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    return(IDASPTFQMR_MEM_FAIL);
  }

  /* Set SPTFQMR parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPTFQMR_MAXL : maxl;
  idasptfqmr_mem->q_maxl = maxl1;

  /* Set default values for the rest of the Sptfqmr parameters */
  idasptfqmr_mem->q_eplifac   = PT05;
  idasptfqmr_mem->q_dqincfac  = ONE;
  idasptfqmr_mem->q_pset      = NULL;
  idasptfqmr_mem->q_psolve    = NULL;
  idasptfqmr_mem->q_pdata     = NULL;
  idasptfqmr_mem->q_jtimes    = IDASptfqmrDQJtimes;
  idasptfqmr_mem->q_jdata     = ida_mem;
  idasptfqmr_mem->q_last_flag = IDASPTFQMR_SUCCESS;

  /* Set setupNonNull to FALSE */
  setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    return(IDASPTFQMR_MEM_FAIL);
  }
  yptemp = N_VClone(vec_tmpl);
  if (yptemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    N_VDestroy(ytemp);
    return(IDASPTFQMR_MEM_FAIL);
  }
  xx = N_VClone(vec_tmpl);
  if (xx == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    return(IDASPTFQMR_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt(N_VDotProd(ytemp, ytemp));

  /* Call SptfqmrMalloc to allocate workspace for Sptfqmr */
  sptfqmr_mem = SptfqmrMalloc(maxl1, vec_tmpl);
  if (sptfqmr_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    N_VDestroy(xx);
    return(IDASPTFQMR_MEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  lmem = idasptfqmr_mem;

  return(IDASPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrSet* and IDASptfqmrGet*
 * -----------------------------------------------------------------
 */

int IDASptfqmrSetMaxl(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  idasptfqmr_mem->q_maxl = (maxl <= 0) ? IDA_SPTFQMR_MAXL : maxl;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrSetEpsLin(void *ida_mem, realtype eplifac)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  if (eplifac < ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_IDAS_NEG_EPLIFAC);
    return(IDASPTFQMR_ILL_INPUT);
  }

  if (eplifac == ZERO)
    idasptfqmr_mem->q_eplifac = PT05;
  else
    idasptfqmr_mem->q_eplifac = eplifac;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrSetIncrementFactor(void *ida_mem, realtype dqincfac)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  if (dqincfac <= ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_IDAS_NEG_DQINCFAC);
    return(IDASPTFQMR_ILL_INPUT);
  }

  idasptfqmr_mem->q_dqincfac = dqincfac;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrSetPreconditioner(void *ida_mem,
				IDASpilsPrecSetupFn pset,
				IDASpilsPrecSolveFn psolve,
				void *prec_data)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  idasptfqmr_mem->q_pset   = pset;
  idasptfqmr_mem->q_psolve = psolve;
  if (psolve != NULL) idasptfqmr_mem->q_pdata = prec_data;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrSetJacTimesVecFn(void *ida_mem, 
			       IDASpilsJacTimesVecFn jtimes,
			       void *jac_data)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  idasptfqmr_mem->q_jtimes = jtimes;
  if (jtimes != NULL) idasptfqmr_mem->q_jdata = jac_data;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS)
{
  IDAMem IDA_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  *lenrwLS = lrw1*13;
  *leniwLS = liw1*13;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumPrecEvals(void *ida_mem, long int *npevals)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *npevals = npe;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumPrecSolves(void *ida_mem, long int *npsolves)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *npsolves = nps;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumLinIters(void *ida_mem, long int *nliters)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *nliters = nli;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *nlcfails = ncfl;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *njvevals = njtimes;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *nrevalsLS = nreSQ;

  return(IDASPTFQMR_SUCCESS);
}

int IDASptfqmrGetLastFlag(void *ida_mem, int *flag)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGTFQMR_SETGET_IDAMEM_NULL);
    return(IDASPTFQMR_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGTFQMR_SETGET_LMEM_NULL);
    return(IDASPTFQMR_LMEM_NULL);
  }
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  *flag = last_flag;

  return(IDASPTFQMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASPTFQMR interface routines
 * -----------------------------------------------------------------
 */

/* Additional readability Replacements */

#define maxl     (idasptfqmr_mem->q_maxl)
#define eplifac  (idasptfqmr_mem->q_eplifac)
#define dqincfac (idasptfqmr_mem->q_dqincfac)
#define psolve   (idasptfqmr_mem->q_psolve)
#define pset     (idasptfqmr_mem->q_pset)
#define pdata    (idasptfqmr_mem->q_pdata)
#define jtimes   (idasptfqmr_mem->q_jtimes)
#define jdata    (idasptfqmr_mem->q_jdata)

static int IDASptfqmrInit(IDAMem IDA_mem)
{
  IDASptfqmrMem idasptfqmr_mem;

  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = 0;
  njtimes = nreSQ = 0;

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = IDASptfqmrDQJtimes;
    jdata = IDA_mem;
  }

  last_flag = IDASPTFQMR_SUCCESS;

  return(0);
}

static int IDASptfqmrSetup(IDAMem IDA_mem, 
			   N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
			   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDASptfqmrMem idasptfqmr_mem;

  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  /* Call user setup routine pset and update counter npe */
  retval = pset(tn, yy_p, yp_p, rr_p, cj, pdata,
                tmp1, tmp2, tmp3);
  npe++;

  last_flag = retval;
  /* Return flag showing success or failure of pset */
  if (retval < 0) return(-1);
  if (retval > 0) return(+1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrSolve
 * -----------------------------------------------------------------
 * Note: The x-scaling and b-scaling arrays are both equal to weight.
 *
 * We set the initial guess, x = 0, then call SptfqmrSolve.
 * We copy the solution x into b, and update the counters nli, nps,
 * and ncfl. If SptfqmrSolve returned nli_inc = 0 (hence x = 0), we
 * take the SPTFQMR vtemp vector (= P_inverse F) as the correction
 * vector instead. Finally, we set the return value according to the
 * success of SptfqmrSolve.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
			   N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDASptfqmrMem idasptfqmr_mem;
  int pretype, nli_inc, nps_inc, retval;
  realtype res_norm;

  idasptfqmr_mem = (IDASptfqmrMem) lmem;


  /* Set SptfqmrSolve convergence test constant epslin, in terms of the
     Newton convergence test constant epsNewt and safety factors. The factor
     sqrt(Neq) assures that the TFQMR convergence test is applied to the
     WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = sqrtN*eplifac*epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  ycur = yy_now;
  ypcur = yp_now;
  rcur = rr_now;

  /* Set SptfqmrSolve inputs pretype and initial guess xx = 0 */  
  pretype = (psolve == NULL) ? PREC_NONE : PREC_LEFT;
  N_VConst(ZERO, xx);
  
  /* Call SptfqmrSolve and copy xx to bb */
  retval = SptfqmrSolve(sptfqmr_mem, IDA_mem, xx, bb, pretype, epslin,
                      IDA_mem, weight, weight, IDASptfqmrAtimes,
                      IDASptfqmrPSolve, &res_norm, &nli_inc, &nps_inc);
  last_flag = retval;
  if (nli_inc == 0) N_VScale(ONE, SPTFQMR_VTEMP(sptfqmr_mem), bb);
  else N_VScale(ONE, xx, bb);
  
  /* Increment counters nli, nps, and return if successful */
  nli += nli_inc;
  nps += nps_inc;

  if (retval == 0) return(0);

  /* If not successful, increment ncfl and return appropriate flag */
  ncfl++;

  if (retval > 0)   return(+1);
  if (retval != -2) return(-1);
  if (resflag > 0)  return(+1);
  return(-1);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrPerf
 * -----------------------------------------------------------------
 * This routine handles performance monitoring specific to the
 * IDASPTFQMR linear solver. When perftask = 0, it saves values of
 * various counters. When perftask = 1, it examines difference
 * quotients in these counters, and depending on their values, it
 * prints up to three warning messages. Messages are printed up to
 * a maximum of 10 times.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrPerf(IDAMem IDA_mem, int perftask)
{
  IDASptfqmrMem idasptfqmr_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  if (perftask == 0) {
    nst0 = nst;  nni0 = nni;  nli0 = nli;
    ncfn0 = ncfn;  ncfl0 = ncfl;  
    nwarn = 0;
    return(0);
  }

  nstd = nst - nst0;  nnid = nni - nni0;
  if (nstd == 0 || nnid == 0) return(0);
  avdim = (realtype) ((nli - nli0)/((realtype) nnid));
  rcfn = (realtype) ((ncfn - ncfn0)/((realtype) nstd));
  rcfl = (realtype) ((ncfl - ncfl0)/((realtype) nnid));
  lavd = (avdim > ((realtype) maxl));
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lavd || lcfn || lcfl)) return(0);
  nwarn++;
  if (nwarn > 10) return(1);
  if (lavd) if (errfp != NULL) fprintf(errfp, MSGTFQMR_AVD_WARN, tn, avdim);
  if (lcfn) if (errfp != NULL) fprintf(errfp, MSGTFQMR_CFN_WARN, tn, rcfn);
  if (lcfl) if (errfp != NULL) fprintf(errfp, MSGTFQMR_CFL_WARN, tn, rcfl);

  return(0);
}

static int IDASptfqmrFree(IDAMem IDA_mem)
{
  IDASptfqmrMem idasptfqmr_mem;

  idasptfqmr_mem = (IDASptfqmrMem) lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(yptemp);
  N_VDestroy(xx);
  SptfqmrFree(sptfqmr_mem);
  free(lmem);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * IDASPTFQMR private functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by calling either the user provided
 * routine or the internal DQ routine.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrAtimes(void *ida_mem, N_Vector v, N_Vector z)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;
  int jtflag;

  IDA_mem = (IDAMem) ida_mem;
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  jtflag = jtimes(tn, ycur, ypcur, rcur, v, z, cj, jdata, ytemp, yptemp);
  njtimes++;

  return(jtflag);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SptfqmrSolve routine
 * and the user's psolve routine. It passes to psolve all required
 * state information from ida_mem. Its return value is the same as
 * that returned by psolve. Note that the generic SPTFQMR solver
 * guarantees that IDASptfqmrPSolve will not be called in the case
 * psolve = NULL.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;
  int retval;

  IDA_mem = (IDAMem) ida_mem;
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  retval = psolve(tn, ycur, ypcur, rcur, r, z, cj, epslin, pdata, ytemp);

  /* This call is counted in nps within the IDASptfqmrSolve routine */

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmrDQJtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by using a difference quotient
 * approximation. The approximation is
 *      Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  where
 *        y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
 *        sigma = sqrt(Neq)*dqincfac.
 * The return value from the call to res is saved in order to set
 * the return flag from IDASptfqmrSolve.
 * -----------------------------------------------------------------
 */

static int IDASptfqmrDQJtimes(realtype tt,
			      N_Vector yy, N_Vector yp, N_Vector rr,
			      N_Vector v, N_Vector Jv, 
			      realtype c_j, void *jac_data, 
			      N_Vector work1, N_Vector work2)
{
  IDAMem IDA_mem;
  IDASptfqmrMem idasptfqmr_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig, siginv;
  int ires;

  /* jac_data is ida_mem */
  IDA_mem = (IDAMem) jac_data;
  idasptfqmr_mem = (IDASptfqmrMem) lmem;

  sig = dqincfac/N_VWrmsNorm(v, ewt);

  /* Rename work1 and work2 for readibility */
  y_tmp  = work1;
  yp_tmp = work2;

  /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
  N_VLinearSum(sig, v, ONE, yy, y_tmp);
  N_VLinearSum(c_j*sig, v, ONE, yp, yp_tmp);

  /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
  ires = res(tt, y_tmp, yp_tmp, Jv, rdata); 
  nreSQ++;
  resflag = ires;
  if (ires != 0) return(ires);

  /* Set Jv to [Jv - rr]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, rr, Jv);

  return(0);
}
