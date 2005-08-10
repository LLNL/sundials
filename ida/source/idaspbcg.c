/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2005-08-10 21:43:22 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the IDA scaled preconditioned
 * Bi-CGSTAB linear solver module, IDASPBCG.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "idaspbcg_impl.h"
#include "sundialsmath.h"

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define PT9  RCONST(0.9)
#define PT05 RCONST(0.05)

#define IDA_SPBCG_MAXL 5

/* IDASPBCG linit, lsetup, lsolve, lperf, and lfree routines */

static int IDASpbcgInit(IDAMem IDA_mem);

static int IDASpbcgSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDASpbcgSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDASpbcgPerf(IDAMem IDA_mem, int perftask);

static int IDASpbcgFree(IDAMem IDA_mem);

/* IDASPBCG Atimes and PSolve routines called by generic SPBCG solver */

static int IDASpbcgAtimes(void *ida_mem, N_Vector v, N_Vector z);

static int IDASpbcgPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

static int IDASpbcgDQJtimes(realtype tt,
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

#define sqrtN     (idaspbcg_mem->b_sqrtN)
#define epslin    (idaspbcg_mem->b_epslin)
#define ytemp     (idaspbcg_mem->b_ytemp)
#define yptemp    (idaspbcg_mem->b_yptemp)
#define xx        (idaspbcg_mem->b_xx)
#define ycur      (idaspbcg_mem->b_ycur)
#define ypcur     (idaspbcg_mem->b_ypcur)
#define rcur      (idaspbcg_mem->b_rcur)
#define resflag   (idaspbcg_mem->b_resflag)
#define npe       (idaspbcg_mem->b_npe)
#define nli       (idaspbcg_mem->b_nli)
#define nps       (idaspbcg_mem->b_nps)
#define ncfl      (idaspbcg_mem->b_ncfl)
#define nst0      (idaspbcg_mem->b_nst0)
#define nni0      (idaspbcg_mem->b_nni0)
#define nli0      (idaspbcg_mem->b_nli0)
#define ncfn0     (idaspbcg_mem->b_ncfn0)
#define ncfl0     (idaspbcg_mem->b_ncfl0)
#define nwarn     (idaspbcg_mem->b_nwarn)
#define njtimes   (idaspbcg_mem->b_njtimes)
#define nreSG     (idaspbcg_mem->b_nreSG)
#define spbcg_mem (idaspbcg_mem->b_spbcg_mem)
#define last_flag (idaspbcg_mem->b_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcg
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDASPBCG linear solver module.
 *
 * IDASpbcg first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDASpbcgInit, IDASpbcgSetup,
 * IDASpbcgSolve, IDASpbcgPerf, and IDASpbcgFree, respectively.
 * It allocates memory for a structure of type IDASpbcgMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem). It then sets the following
 * fields in the IDASpbcgMemRec structure:
 *   b_maxl      = IDA_SPBCG_MAXL  if maxl <= 0
 *               = maxl            if maxl > 0
 *   b_eplifac   = PT05
 *   b_dqincfac  = ONE
 *   b_pdata     = NULL
 *   b_pset      = NULL
 *   b_psolve    = NULL
 *   b_jtimes    = IDASpbcgDQJtimes
 *   b_jdata     = ida_mem
 *   b_last_flag = IDASPBCG_SUCCESS
 * Finally, IDASpbcg allocates memory for ytemp, yptemp, and xx, and
 * calls SpbcgMalloc to allocate memory for the Spbcg solver.
 *
 * The return value of IDASpbcg is:
 *   IDASPBCG_SUCCESS   =  0 if successful
 *   IDASPBCG_MEM_FAIL  = -1 if IDA_mem is NULL or a memory
 *                           allocation failed
 *   IDASPBCG_ILL_INPUT = -2 if a required vector operation is not
 *                           implemented.
 * -----------------------------------------------------------------
 */

int IDASpbcg(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;
  int flag, maxl1;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_BAD_NVECTOR);
    return(IDASPBCG_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  linit  = IDASpbcgInit;
  lsetup = IDASpbcgSetup;
  lsolve = IDASpbcgSolve;
  lperf  = IDASpbcgPerf;
  lfree  = IDASpbcgFree;

  /* Get memory for IDASpbcgMemRec */
  idaspbcg_mem = (IDASpbcgMem) malloc(sizeof(IDASpbcgMemRec));
  if (idaspbcg_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    return(IDASPBCG_MEM_FAIL);
  }

  /* Set SPBCG parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPBCG_MAXL : maxl;
  idaspbcg_mem->b_maxl = maxl1;

  /* Set default values for the rest of the Spbcg parameters */
  idaspbcg_mem->b_eplifac   = PT05;
  idaspbcg_mem->b_dqincfac  = ONE;
  idaspbcg_mem->b_pset      = NULL;
  idaspbcg_mem->b_psolve    = NULL;
  idaspbcg_mem->b_pdata     = NULL;
  idaspbcg_mem->b_jtimes    = IDASpbcgDQJtimes;
  idaspbcg_mem->b_jdata     = ida_mem;
  idaspbcg_mem->b_last_flag = IDASPBCG_SUCCESS;

  /* Set setupNonNull to FALSE */
  setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    return(IDASPBCG_MEM_FAIL);
  }
  yptemp = N_VClone(vec_tmpl);
  if (yptemp == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    N_VDestroy(ytemp);
    return(IDASPBCG_MEM_FAIL);
  }
  xx = N_VClone(vec_tmpl);
  if (xx == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    return(IDASPBCG_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt(N_VDotProd(ytemp, ytemp));

  /* Call SpbcgMalloc to allocate workspace for Spbcg */
  spbcg_mem = SpbcgMalloc(maxl1, vec_tmpl);
  if (spbcg_mem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    N_VDestroy(xx);
    return(IDASPBCG_MEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  lmem = idaspbcg_mem;

  return(IDASPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgSet* and IDASpbcgGet*
 * -----------------------------------------------------------------
 */

int IDASpbcgSetEpsLin(void *ida_mem, realtype eplifac)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  if (eplifac < ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_IDAS_NEG_EPLIFAC);
    return(IDASPBCG_ILL_INPUT);
  }

  if (eplifac == ZERO)
    idaspbcg_mem->b_eplifac = PT05;
  else
    idaspbcg_mem->b_eplifac = eplifac;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgSetIncrementFactor(void *ida_mem, realtype dqincfac)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  if (dqincfac <= ZERO) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_IDAS_NEG_DQINCFAC);
    return(IDASPBCG_ILL_INPUT);
  }

  idaspbcg_mem->b_dqincfac = dqincfac;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
                              IDASpilsPrecSolveFn psolve, void *prec_data)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  idaspbcg_mem->b_pset = pset;
  idaspbcg_mem->b_psolve = psolve;
  if (psolve != NULL) idaspbcg_mem->b_pdata = prec_data;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgSetJacTimesVecFn(void *ida_mem, 
                             IDASpilsJacTimesVecFn jtimes, void *jac_data)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  idaspbcg_mem->b_jtimes = jtimes;
  if (jtimes != NULL) idaspbcg_mem->b_jdata = jac_data;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetWorkSpace(void *ida_mem, long int *lenrwSG, long int *leniwSG)
{
  IDAMem IDA_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }

  *lenrwSG = lrw1 * 10;
  *leniwSG = liw1 * 10;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumPrecEvals(void *ida_mem, long int *npevals)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *npevals = npe;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumPrecSolves(void *ida_mem, long int *npsolves)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *npsolves = nps;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumLinIters(void *ida_mem, long int *nliters)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *nliters = nli;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *nlcfails = ncfl;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *njvevals = njtimes;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetNumResEvals(void *ida_mem, long int *nrevalsSG)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *nrevalsSG = nreSG;

  return(IDASPBCG_SUCCESS);
}

int IDASpbcgGetLastFlag(void *ida_mem, int *flag)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSGBCG_SETGET_IDAMEM_NULL);
    return(IDASPBCG_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if (errfp != NULL) fprintf(errfp, MSGBCG_SETGET_LMEM_NULL);
    return(IDASPBCG_LMEM_NULL);
  }
  idaspbcg_mem = (IDASpbcgMem) lmem;

  *flag = last_flag;

  return(IDASPBCG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASPBCG interface routines
 * -----------------------------------------------------------------
 */

/* Additional readability Replacements */

#define maxl     (idaspbcg_mem->b_maxl)
#define eplifac  (idaspbcg_mem->b_eplifac)
#define dqincfac (idaspbcg_mem->b_dqincfac)
#define psolve   (idaspbcg_mem->b_psolve)
#define pset     (idaspbcg_mem->b_pset)
#define pdata    (idaspbcg_mem->b_pdata)
#define jtimes   (idaspbcg_mem->b_jtimes)
#define jdata    (idaspbcg_mem->b_jdata)

static int IDASpbcgInit(IDAMem IDA_mem)
{
  IDASpbcgMem idaspbcg_mem;

  idaspbcg_mem = (IDASpbcgMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = 0;
  njtimes = nreSG = 0;

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = IDASpbcgDQJtimes;
    jdata = IDA_mem;
  }

  last_flag = IDASPBCG_SUCCESS;

  return(0);
}

static int IDASpbcgSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDASpbcgMem idaspbcg_mem;

  idaspbcg_mem = (IDASpbcgMem) lmem;

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
 * Function : IDASpbcgSolve
 * -----------------------------------------------------------------
 * Note: The x-scaling and b-scaling arrays are both equal to weight.
 *
 * We set the initial guess, x = 0, then call SpbcgSolve.
 * We copy the solution x into b, and update the counters nli, nps,
 * and ncfl. If SpbcgSolve returned nli_inc = 0 (hence x = 0), we
 * take the SPBCG vtemp vector (= P_inverse F) as the correction
 * vector instead. Finally, we set the return value according to the
 * success of SpbcgSolve.
 * -----------------------------------------------------------------
 */

static int IDASpbcgSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDASpbcgMem idaspbcg_mem;
  int pretype, nli_inc, nps_inc, retval;
  realtype res_norm;

  idaspbcg_mem = (IDASpbcgMem) lmem;


  /* Set SpbcgSolve convergence test constant epslin, in terms of the
     Newton convergence test constant epsNewt and safety factors. The factor
     sqrt(Neq) assures that the Bi-CGSTAB convergence test is applied to the
     WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = sqrtN*eplifac*epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  ycur = yy_now;
  ypcur = yp_now;
  rcur = rr_now;

  /* Set SpbcgSolve inputs pretype and initial guess xx = 0 */  
  pretype = (psolve == NULL) ? PREC_NONE : PREC_LEFT;
  N_VConst(ZERO, xx);
  
  /* Call SpbcgSolve and copy xx to bb */
  retval = SpbcgSolve(spbcg_mem, IDA_mem, xx, bb, pretype, epslin,
                      IDA_mem, weight, weight, IDASpbcgAtimes,
                      IDASpbcgPSolve, &res_norm, &nli_inc, &nps_inc);
  last_flag = retval;
  if (nli_inc == 0) N_VScale(ONE, SPBCG_VTEMP(spbcg_mem), bb);
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
 * Function : IDASpbcgPerf
 * -----------------------------------------------------------------
 * This routine handles performance monitoring specific to the
 * IDASPBCG linear solver. When perftask = 0, it saves values of
 * various counters. When perftask = 1, it examines difference
 * quotients in these counters, and depending on their values, it
 * prints up to three warning messages. Messages are printed up to
 * a maximum of 10 times.
 * -----------------------------------------------------------------
 */

static int IDASpbcgPerf(IDAMem IDA_mem, int perftask)
{
  IDASpbcgMem idaspbcg_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idaspbcg_mem = (IDASpbcgMem) lmem;

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
  if (lavd) if (errfp != NULL) fprintf(errfp, MSGBCG_AVD_WARN, tn, avdim);
  if (lcfn) if (errfp != NULL) fprintf(errfp, MSGBCG_CFN_WARN, tn, rcfn);
  if (lcfl) if (errfp != NULL) fprintf(errfp, MSGBCG_CFL_WARN, tn, rcfl);

  return(0);
}

static int IDASpbcgFree(IDAMem IDA_mem)
{
  IDASpbcgMem idaspbcg_mem;

  idaspbcg_mem = (IDASpbcgMem) lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(xx);
  SpbcgFree(spbcg_mem);
  free(lmem);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * IDASPBCG private functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by calling either the user provided
 * routine or the internal DQ routine.
 * -----------------------------------------------------------------
 */

static int IDASpbcgAtimes(void *ida_mem, N_Vector v, N_Vector z)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;
  int jtflag;

  IDA_mem = (IDAMem) ida_mem;
  idaspbcg_mem = (IDASpbcgMem) lmem;

  jtflag = jtimes(tn, ycur, ypcur, rcur, v, z, cj, jdata, ytemp, yptemp);
  njtimes++;

  return(jtflag);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SpbcgSolve routine
 * and the user's psolve routine. It passes to psolve all required
 * state information from ida_mem. Its return value is the same as
 * that returned by psolve. Note that the generic SPBCG solver
 * guarantees that IDASpbcgPSolve will not be called in the case
 * psolve = NULL.
 * -----------------------------------------------------------------
 */

static int IDASpbcgPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;
  int retval;

  IDA_mem = (IDAMem) ida_mem;
  idaspbcg_mem = (IDASpbcgMem) lmem;

  retval = psolve(tn, ycur, ypcur, rcur, r, z, cj, epslin, pdata, ytemp);

  /* This call is counted in nps within the IDASpbcgSolve routine */

  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgDQJtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Jv, where
 * J is the system Jacobian, by using a difference quotient
 * approximation. The approximation is
 *      Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  where
 *        y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
 *        sigma = sqrt(Neq)*dqincfac.
 * The return value from the call to res is saved in order to set
 * the return flag from IDASpbcgSolve.
 * -----------------------------------------------------------------
 */

static int IDASpbcgDQJtimes(realtype tt,
                            N_Vector yy, N_Vector yp, N_Vector rr,
                            N_Vector v, N_Vector Jv, 
                            realtype c_j, void *jac_data, 
                            N_Vector work1, N_Vector work2)
{
  IDAMem IDA_mem;
  IDASpbcgMem idaspbcg_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig, siginv;
  int ires;

  /* jac_data is ida_mem */
  IDA_mem = (IDAMem) jac_data;
  idaspbcg_mem = (IDASpbcgMem) lmem;

  sig = dqincfac/N_VWrmsNorm(v, ewt);

  /* Rename work1 and work2 for readibility */
  y_tmp  = work1;
  yp_tmp = work2;

  /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
  N_VLinearSum(sig, v, ONE, yy, y_tmp);
  N_VLinearSum(c_j*sig, v, ONE, yp, yp_tmp);

  /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
  ires = res(tt, y_tmp, yp_tmp, Jv, rdata); 
  nreSG++;
  resflag = ires;
  if (ires != 0) return(ires);

  /* Set Jv to [Jv - rr]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, rr, Jv);

  return(0);
}
