/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2004-07-22 23:01:09 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for the IDA Scaled              
 * Preconditioned GMRES linear solver module, IDASPGMR.            
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "ida_impl.h"
#include "idaspgmr_impl.h"
#include "spgmr.h"
#include "iterative.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"

/* Error Messages */

#define IDASPGMR               "IDASpgmr-- "

#define MSG_IDAMEM_NULL        IDASPGMR "Integrator memory is NULL.\n\n"

#define MSG_MEM_FAIL           IDASPGMR "A memory request failed.\n\n"

#define MSG_BAD_NVECTOR        IDASPGMR "A required vector operation is not implemented.\n\n"

#define MSG_SETGET_IDAMEM_NULL "IDASpgmrSet*/IDASpgmrGet*-- Integrator memory is NULL. \n\n"

#define MSG_SETGET_LMEM_NULL   "IDASpgmrSet*/IDASpgmrGet*-- IDASPGMR memory is NULL. \n\n"

#define MSG_BAD_GSTYPE_1       "IDASpgmrSetGSType-- gstype = %d illegal.\n"
#define MSG_BAD_GSTYPE_2       "The legal values are MODIFIED_GS = %d and "
#define MSG_BAD_GSTYPE_3       "CLASSICAL_GS = %d.\n\n"
#define MSG_BAD_GSTYPE         MSG_BAD_GSTYPE_1 MSG_BAD_GSTYPE_2 MSG_BAD_GSTYPE_3

#define MSG_IDAS_NEG_MAXRS     "IDASpgmrSetMaxRestarts-- maxrs<0 illegal. \n\n"

#define MSG_IDAS_NEG_EPLIFAC   "IDASpgmrSetEpsLin-- eplifac < 0.0 illegal. \n\n"

#define MSG_IDAS_NEG_DQINCFAC  "IDASpgmrSetIncrementFactor-- dqincfac < 0.0 illegal. \n\n"

/* Warning Messages */

#define IDASpgmr_PERF  "IDASpgmrPerf--"

#define MSG_WARN1      IDASpgmr_PERF "Warning. Poor iterative algorithm"
#define MSG_WARN2      " performance at t = %e\n"

#define MSG_AVD_WARN3  "Average number of linear iterations is %e.\n\n"
#define MSG_AVD_WARN   MSG_WARN1 MSG_WARN2 MSG_AVD_WARN3

#define MSG_CFN_WARN3  "Nonlinear convergence failure rate is %e.\n\n"
#define MSG_CFN_WARN   MSG_WARN1 MSG_WARN2 MSG_CFN_WARN3

#define MSG_CFL_WARN3  "Linear convergence failure rate is %e.\n\n"
#define MSG_CFL_WARN   MSG_WARN1 MSG_WARN2 MSG_CFL_WARN3


/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define PT9          RCONST(0.9)
#define PT05         RCONST(0.05)

#define IDA_SPGMR_MAXL    5
#define IDA_SPGMR_MAXRS   5


/* IDASPGMR linit, lsetup, lsolve, lperf, and lfree routines */

static int IDASpgmrInit(IDAMem IDA_mem);

static int IDASpgmrSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1,
                         N_Vector tempv2, N_Vector tempv3);

static int IDASpgmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector ynow, N_Vector ypnow, N_Vector rnow);

static int IDASpgmrPerf(IDAMem IDA_mem, int perftask);

static int IDASpgmrFree(IDAMem IDA_mem);


/* IDASPGMR Atimes and PSolve routines called by generic SPGMR solver */

static int IDASpgmrAtimes(void *ida_mem, N_Vector v, N_Vector z);

static int IDASpgmrPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

static int IDASpgmrDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                            N_Vector yy, N_Vector yp, N_Vector rr,
                            realtype c_j, void *jac_data, 
                            N_Vector work1, N_Vector work2);
/* Readability Replacements */

#define lrw1    (IDA_mem->ida_lrw1)
#define liw1    (IDA_mem->ida_liw1)
#define nst     (IDA_mem->ida_nst)
#define tn      (IDA_mem->ida_tn)
#define cj      (IDA_mem->ida_cj)
#define epsNewt (IDA_mem->ida_epsNewt)
#define nre     (IDA_mem->ida_nre)
#define res     (IDA_mem->ida_res)
#define rdata   (IDA_mem->ida_rdata)
#define ewt     (IDA_mem->ida_ewt)
#define errfp   (IDA_mem->ida_errfp)
#define iopt    (IDA_mem->ida_iopt)
#define linit   (IDA_mem->ida_linit)
#define lsetup  (IDA_mem->ida_lsetup)
#define lsolve  (IDA_mem->ida_lsolve)
#define lperf   (IDA_mem->ida_lperf)
#define lfree   (IDA_mem->ida_lfree)
#define lmem    (IDA_mem->ida_lmem)
#define nni     (IDA_mem->ida_nni)
#define ncfn    (IDA_mem->ida_ncfn)
#define setupNonNull (IDA_mem->ida_setupNonNull)
#define vec_tmpl     (IDA_mem->ida_tempv1)

#define sqrtN   (idaspgmr_mem->g_sqrtN)
#define epslin  (idaspgmr_mem->g_epslin)
#define ytemp   (idaspgmr_mem->g_ytemp)
#define yptemp  (idaspgmr_mem->g_yptemp)
#define xx      (idaspgmr_mem->g_xx)
#define ycur    (idaspgmr_mem->g_ycur)
#define ypcur   (idaspgmr_mem->g_ypcur)
#define rcur    (idaspgmr_mem->g_rcur)
#define resflag (idaspgmr_mem->g_resflag)
#define npe     (idaspgmr_mem->g_npe)
#define nli     (idaspgmr_mem->g_nli)
#define nps     (idaspgmr_mem->g_nps)
#define ncfl    (idaspgmr_mem->g_ncfl)
#define nst0    (idaspgmr_mem->g_nst0)
#define nni0    (idaspgmr_mem->g_nni0)
#define nli0    (idaspgmr_mem->g_nli0)
#define ncfn0   (idaspgmr_mem->g_ncfn0)
#define ncfl0   (idaspgmr_mem->g_ncfl0)
#define nwarn   (idaspgmr_mem->g_nwarn)
#define njtimes (idaspgmr_mem->g_njtimes)
#define nreSG   (idaspgmr_mem->g_nreSG)

#define spgmr_mem (idaspgmr_mem->g_spgmr_mem)


/*************** IDASpgmr *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the IDASPGMR linear solver module.  

 IDASpgmr first calls the existing lfree routine if this is not NULL.
 It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 ida_lfree fields in (*IDA_mem) to be IDASpgmrInit, IDASpgmrSetup,
 IDASpgmrSolve, IDASpgmrPerf, and IDASpgmrFree, respectively.
 It allocates memory for a structure of type IDASpgmrMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.
 It sets setupNonNull in (*IDA_mem).  It then sets the following
 fields in the IDASpgmrMemRec structure:
   g_gstype   = gstype
   g_maxl     = MIN(Neq,IDA_SPGMR_MAXL) if maxl <= 0,  else MIN(Neq,maxl)
   g_maxrs    = 0 if maxrs < 0,  MIN(5,Neq/g_maxl) if maxrs = 0, and
                MIN(maxrs,Neq/g_maxl) if maxrs > 0.
   g_eplifac  = 0.05 if eplifac = 0.0,  else eplifac
   g_dqincfac = 1.0 if dqincfac = 0.0,  else dqincfac
   g_pdata    = NULL
   g_pset     = NULL
   g_psolve   = NULL
   g_jtimes   = NULL
   g_jdata    = NULL
 Finally, IDASpgmr allocates memory for ytemp, yptemp, and xx, and
 calls SpgmrMalloc to allocate memory for the Spgmr solver.

 The return value of IDASpgmr is:
   SUCCESS       = 0  if successful
   LMEM_FAIL     = -1 if IDA_mem is NULL or a memory allocation failed
   LIN_ILL_INPUT = -2 if the gstype argument is illegal.

**********************************************************************/

int IDASpgmr(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  int flag, maxl1;
  booleantype nvectorOK;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_IDAMEM_NULL);
    return(LMEM_FAIL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if(vec_tmpl->ops->nvdotprod == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(LIN_ILL_INPUT);
  }

  if (lfree != NULL) flag = lfree(ida_mem);

  /* Set five main function fields in ida_mem */
  linit  = IDASpgmrInit;
  lsetup = IDASpgmrSetup;
  lsolve = IDASpgmrSolve;
  lperf  = IDASpgmrPerf;
  lfree  = IDASpgmrFree;

  /* Get memory for IDASpgmrMemRec */
  idaspgmr_mem = (IDASpgmrMem) malloc(sizeof(IDASpgmrMemRec));
  if (idaspgmr_mem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }

  /* Set SPGMR parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPGMR_MAXL : maxl;
  idaspgmr_mem->g_maxl     = maxl1;

  /* Set default values for the rest of the Spgmr parameters */
  idaspgmr_mem->g_gstype   = MODIFIED_GS;
  idaspgmr_mem->g_maxrs    = IDA_SPGMR_MAXRS;
  idaspgmr_mem->g_eplifac  = PT05;
  idaspgmr_mem->g_dqincfac = ONE;
  idaspgmr_mem->g_pset     = NULL;
  idaspgmr_mem->g_psolve   = NULL;
  idaspgmr_mem->g_pdata    = NULL;
  idaspgmr_mem->g_jtimes   = IDASpgmrDQJtimes;
  idaspgmr_mem->g_jdata    = ida_mem;

  /* Set setupNonNull to FALSE */
  setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */
  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }
  yptemp = N_VClone(vec_tmpl);
  if (yptemp == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    N_VDestroy(ytemp);
    return(LMEM_FAIL);
  }
  xx = N_VClone(vec_tmpl);
  if (xx == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    return(LMEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt( N_VDotProd(ytemp, ytemp) );

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(maxl1, vec_tmpl);
  if (spgmr_mem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(yptemp);
    N_VDestroy(xx);
    return(LMEM_FAIL);
  }

  /* Attach linear solver memory to the integrator memory */
  lmem = idaspgmr_mem;

  return(SUCCESS);
}

/*************** IDASpgmrSetGSType ************************************/

int IDASpgmrSetGSType(void *ida_mem, int gstype)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
    return(LIN_ILL_INPUT);
  }

  idaspgmr_mem->g_gstype = gstype;

  return(SUCCESS);
}

/*************** IDASpgmrSetMaxRestarts ******************************/

int IDASpgmrSetMaxRestarts(void *ida_mem, int maxrs)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Check for legal maxrs */
  if (maxrs < 0) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_MAXRS);
    return(LIN_ILL_INPUT);
  }

  idaspgmr_mem->g_maxrs = maxrs;

  return(SUCCESS);
}

/*************** IDASpgmrSetEpsLin ***********************************/

int IDASpgmrSetEpsLin(void *ida_mem, realtype eplifac)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Check for legal maxrs */
  if (eplifac < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_EPLIFAC);
    return(LIN_ILL_INPUT);
  }

  idaspgmr_mem->g_eplifac = eplifac;

  return(SUCCESS);
}

/*************** IDASpgmrSetIncrementFactor **************************/

int IDASpgmrSetIncrementFactor(void *ida_mem, realtype dqincfac)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Check for legal maxrs */
  if (dqincfac < ZERO) {
    if(errfp!=NULL) fprintf(errfp, MSG_IDAS_NEG_DQINCFAC);
    return(LIN_ILL_INPUT);
  }

  idaspgmr_mem->g_dqincfac = dqincfac;

  return(SUCCESS);
}

/*************** IDASpgmrSetPrecSetupFn *******************************/

int IDASpgmrSetPrecSetupFn(void *ida_mem, IDASpgmrPrecSetupFn pset)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  idaspgmr_mem->g_pset = pset;

  return(SUCCESS);
}

/*************** IDASpgmrSetPrecSolveFn *******************************/

int IDASpgmrSetPrecSolveFn(void *ida_mem, IDASpgmrPrecSolveFn psolve)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  idaspgmr_mem->g_psolve = psolve;

  return(SUCCESS);
}

/*************** IDASpgmrSetPrecData **********************************/

int IDASpgmrSetPrecData(void *ida_mem, void *pdata)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  idaspgmr_mem->g_pdata = pdata;

  return(SUCCESS);
}

/*************** IDASpgmrSetJacTimesVecFn *****************************/

int IDASpgmrSetJacTimesVecFn(void *ida_mem, IDASpgmrJacTimesVecFn jtimes)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  idaspgmr_mem->g_jtimes = jtimes;

  return(SUCCESS);
}

/*************** IDASpgmrSetJacData ***********************************/

int IDASpgmrSetJacData(void *ida_mem, void *jdata)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  idaspgmr_mem->g_jdata = jdata;

  return(SUCCESS);
}

/*************** IDASpgmrGetIntWorkSpace ******************************/

int IDASpgmrGetIntWorkSpace(void *ida_mem, long int *leniwSG)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  int maxl;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  maxl = idaspgmr_mem->g_maxl;
  *leniwSG = liw1*(maxl + 6);

  return(OKAY);
}

/*************** IDASpgmrGetRealWorkSpace *****************************/

int IDASpgmrGetRealWorkSpace(void *ida_mem, long int *lenrwSG)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  int maxl;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  maxl = idaspgmr_mem->g_maxl;
  *lenrwSG = lrw1*(maxl + 6) + maxl*(maxl + 4) + 1;

  return(OKAY);
}

/*************** IDASpgmrGetNumPrecEvals ******************************/

int IDASpgmrGetNumPrecEvals(void *ida_mem, long int *npevals)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *npevals = npe;

  return(OKAY);
}

/*************** IDASpgmrGetNumPrecSolves *****************************/

int IDASpgmrGetNumPrecSolves(void *ida_mem, long int *npsolves)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *npsolves = nps;

  return(OKAY);
}

/*************** IDASpgmrGetNumLinIters *******************************/

int IDASpgmrGetNumLinIters(void *ida_mem, long int *nliters)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *nliters = nli;

  return(OKAY);
}

/*************** IDASpgmrGetNumConvFails ******************************/

int IDASpgmrGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *nlcfails = ncfl;

  return(OKAY);
}

/*************** IDASpgmrGetNumJtimesEvals ****************************/

int IDASpgmrGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *njvevals = njtimes;

  return(OKAY);
}

/*************** IDASpgmrGetNumResEvals *******************************/

int IDASpgmrGetNumResEvals(void *ida_mem, long int *nrevalsSG)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    fprintf(stderr, MSG_SETGET_IDAMEM_NULL);
    return(LIN_NO_MEM);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  idaspgmr_mem = (IDASpgmrMem) lmem;

  *nrevalsSG = nreSG;

  return(OKAY);
}


/* Additional readability Replacements */

#define gstype   (idaspgmr_mem->g_gstype)
#define maxl     (idaspgmr_mem->g_maxl)
#define maxrs    (idaspgmr_mem->g_maxrs)
#define eplifac  (idaspgmr_mem->g_eplifac)
#define dqincfac (idaspgmr_mem->g_dqincfac)
#define psolve   (idaspgmr_mem->g_psolve)
#define pset     (idaspgmr_mem->g_pset)
#define pdata    (idaspgmr_mem->g_pdata)
#define jtimes   (idaspgmr_mem->g_jtimes)
#define jdata    (idaspgmr_mem->g_jdata)

/*************** IDASpgmrInit *****************************************

 This routine does remaining initializations specific to the IDASPGMR 
 linear solver.  It returns LINIT_OK = 0.

**********************************************************************/

static int IDASpgmrInit(IDAMem IDA_mem)
{
  IDASpgmrMem idaspgmr_mem;

  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = 0;
  njtimes = nreSG = 0;

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = IDASpgmrDQJtimes;
    jdata = IDA_mem;
  }

  return(LINIT_OK);
}

/*************** IDASpgmrSetup ****************************************

 This routine calls the user's preconditioner setup routine pset,
 and updates the counter npe.
 The return value is either
     SUCCESS = 0           if successful,
     LSETUP_ERROR_RECVR    if pset failed recoverably, or
     LSETUP_ERROR_NONRECVR if pset failed unrecoverably.

**********************************************************************/

static int IDASpgmrSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1, N_Vector tempv2,
                         N_Vector tempv3)
{
  int retval;
  IDASpgmrMem idaspgmr_mem;

  idaspgmr_mem = (IDASpgmrMem) lmem;

  /* Call user setup routine pset and update counter npe. */
  retval = pset(tn, yyp, ypp, resp, cj, pdata,
                tempv1, tempv2, tempv3);
  npe++;

  /* Return flag showing success or failure of pset. */
  if (retval < 0) return(LSETUP_ERROR_NONRECVR);
  if (retval > 0) return(LSETUP_ERROR_RECVR);
  return(SUCCESS);
}

/*************** IDASpgmrSolve ****************************************

 This routine handles the call to the generic SPGMR solver SpgmrSolve
 for the solution of the linear system Ax = b.

 The x-scaling and b-scaling arrays are both equal to weight.

 We set the initial guess, x = 0, then call SpgmrSolve.  
 We copy the solution x into b, and update the counters nli, nps, ncfl.
 If SpgmrSolve returned nli_inc = 0 (hence x = 0), we take the SPGMR
 vtemp vector (= P_inverse F) as the correction vector instead.
 Finally, we set the return value according to the success of SpgmrSolve.

**********************************************************************/

static int IDASpgmrSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector ynow, N_Vector ypnow, N_Vector rnow)
{
  IDASpgmrMem idaspgmr_mem;
  int pretype, nli_inc, nps_inc, retval;
  realtype res_norm;
  
  idaspgmr_mem = (IDASpgmrMem) lmem;


  /* Set SpgmrSolve convergence test constant epslin, in terms of the
    Newton convergence test constant epsNewt and safety factors.  The factor 
    sqrt(Neq) assures that the GMRES convergence test is applied to the
    WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = sqrtN*eplifac*epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  ycur = ynow;
  ypcur = ypnow;
  rcur = rnow;

  /* Set SpgmrSolve inputs pretype and initial guess xx = 0. */  
  pretype = (psolve == NULL) ? NONE : LEFT;
  N_VConst(ZERO, xx);
  
  /* Call SpgmrSolve and copy xx to bb. */
  retval = SpgmrSolve(spgmr_mem, IDA_mem, xx, bb, pretype, gstype, epslin,
                      maxrs, IDA_mem, weight, weight, IDASpgmrAtimes,
                      IDASpgmrPSolve, &res_norm, &nli_inc, &nps_inc);

  if (nli_inc == 0) N_VScale(ONE, SPGMR_VTEMP(spgmr_mem), bb);
  else N_VScale(ONE, xx, bb);
  
  /* Increment counters nli, nps, and return if successful. */
  nli += nli_inc;
  nps += nps_inc;

  if (retval == SUCCESS) return(SUCCESS);

  /* If not successful, increment ncfl and return appropriate flag. */
  ncfl++;

  if (retval > 0)   return(LSOLVE_ERROR_RECVR);
  if (retval != -2) return(LSOLVE_ERROR_NONRECVR);
  if (resflag > 0)  return(LSOLVE_ERROR_RECVR);
  return(LSOLVE_ERROR_NONRECVR);

}

/*************** IDASpgmrPerf *****************************************

 This routine handles performance monitoring specific to the IDASPGMR
 linear solver.  When perftask = 0, it saves values of various counters.
 When perftask = 1, it examines difference quotients in these counters,
 and depending on their values, it prints up to three warning messages.
 Messages are printed up to a maximum of 10 times.

**********************************************************************/

static int IDASpgmrPerf(IDAMem IDA_mem, int perftask)
{
  IDASpgmrMem idaspgmr_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idaspgmr_mem = (IDASpgmrMem) lmem;

  if (perftask == 0) {
    nst0 = nst;  nni0 = nni;  nli0 = nli;
    ncfn0 = ncfn;  ncfl0 = ncfl;  
    nwarn = 0;
    return 0;
  }

  nstd = nst - nst0;  nnid = nni - nni0;
  if (nstd == 0 || nnid == 0) return 0;
  avdim = (nli - nli0)/( (realtype) nnid);
  rcfn = (ncfn - ncfn0)/( (realtype) nstd);
  rcfl = (ncfl - ncfl0)/( (realtype) nnid);
  lavd = (avdim > ( (realtype) maxl ) );
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lavd || lcfn || lcfl)) return 0;
  nwarn++;
  if (nwarn > 10) return 1;
  if (lavd) if(errfp!=NULL) fprintf(errfp, MSG_AVD_WARN, tn, avdim);
  if (lcfn) if(errfp!=NULL) fprintf(errfp, MSG_CFN_WARN, tn, rcfn);
  if (lcfl) if(errfp!=NULL) fprintf(errfp, MSG_CFL_WARN, tn, rcfl);

  return 0;
}

/*************** IDASpgmrFree *****************************************

 This routine frees memory specific to the IDASPGMR linear solver.

**********************************************************************/

static int IDASpgmrFree(IDAMem IDA_mem)
{
  IDASpgmrMem idaspgmr_mem;

  idaspgmr_mem = (IDASpgmrMem) lmem;
  
  N_VDestroy(ytemp);
  N_VDestroy(xx);
  SpgmrFree(spgmr_mem);
  free(lmem);

  return 0;
}

/***************** IDASpgmrAtimes *************************************

 This routine generates the matrix-vector product z = Jv, where
 J is the system Jacobian, by calling either the user provided
 routine or the internal DQ routine.

**********************************************************************/

static int IDASpgmrAtimes(void *ida_mem, N_Vector v, N_Vector z)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  int jtflag;

  IDA_mem = (IDAMem) ida_mem;
  idaspgmr_mem = (IDASpgmrMem) lmem;

  jtflag = jtimes(v, z, tn, ycur, ypcur, rcur, cj, jdata, ytemp, yptemp);
  njtimes++;

  return(jtflag);
}

/*************** IDASpgmrPSolve ***************************************

 This routine interfaces between the generic SpgmrSolve routine and
 the user's psolve routine.  It passes to psolve all required state 
 information from ida_mem.  Its return value is the same as that
 returned by psolve.  Note that the generic SPGMR solver guarantees
 that IDASpgmrPSolve will not be called in the case psolve = NULL.

**********************************************************************/

static int IDASpgmrPSolve(void *ida_mem, N_Vector r, N_Vector z, int lr)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  int retval;

  IDA_mem = (IDAMem) ida_mem;
  idaspgmr_mem = (IDASpgmrMem) lmem;

  retval = psolve(tn, ycur, ypcur, rcur, r, z, cj, epslin, pdata, ytemp);

  /* This call is counted in nps within the IDASpgmrSolve routine */

  return(retval);

}

/******************* IDASpgmrDQJtimes *********************************

 This routine generates the matrix-vector product z = Jv, where
 J is the system Jacobian, by using a difference quotient approximation.
 The approximation is 
      Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  where
        y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
        sigma = sqrt(Neq)*dqincfac.
 The return value from the call to res is saved in order to set the
 return flag from IDASpgmrSolve.

**********************************************************************/

static int IDASpgmrDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                            N_Vector yy, N_Vector yp, N_Vector rr,
                            realtype c_j, void *jac_data, 
                            N_Vector work1, N_Vector work2)
{
  IDAMem IDA_mem;
  IDASpgmrMem idaspgmr_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig, siginv;
  int ires;

  /* jac_data is ida_mem */
  IDA_mem = (IDAMem) jac_data;
  idaspgmr_mem = (IDASpgmrMem) lmem;

  sig = sqrtN*dqincfac;

  /* Rename work1 and work2 for readibility */
  y_tmp  = work1;
  yp_tmp = work2;

  /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
  N_VLinearSum(sig, v, ONE, yy, ytemp);
  N_VLinearSum(c_j*sig, v, ONE, yp, yptemp);

  /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
  ires = res(t, y_tmp, yp_tmp, Jv, rdata); 
  nreSG++;
  resflag = ires;
  if (ires != 0) return(ires);

  /* Set Jv to [Jv - rr]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, rr, Jv);

  return(0);

}
