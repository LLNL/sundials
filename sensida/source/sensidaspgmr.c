/******************************************************************
 *                                                                *
 * File          : sensidaspgmr.c                                 *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 21 March 2001                                  *
 *----------------------------------------------------------------*
 * This is the implementation file for the IDA Scaled             *
 * Preconditioned GMRES linear solver module, SensIDASPGMR.       *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "sensidaspgmr.h"
#include "ida.h"
#include "spgmr.h"
#include "iterativ.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "sensida.h"

/* Error Messages */

#define MSG_BAD_GSTYPE_1  "SensIDASpgmr-- gstype = %d illegal.\n"
#define MSG_BAD_GSTYPE_2  "The legal values are MODIFIED_GS = %d and "
#define MSG_BAD_GSTYPE_3  "CLASSICAL_GS = %d.\n\n"
#define MSG_BAD_GSTYPE     MSG_BAD_GSTYPE_1 MSG_BAD_GSTYPE_2 MSG_BAD_GSTYPE_3

#define MSG_MEM_FAIL      "SensIDASpgmrInit-- A memory request failed.\n\n"


/* Warning Messages */

#define IDASpgmr_PERF  "SensIDASpgmrPerf--"

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


/******************************************************************
 *                                                                *           
 * Types : IDASpgmrMemRec, IDASpgmrMem                            *
 *----------------------------------------------------------------*
 * The type IDASpgmrMem is pointer to an IDASpgmrMemRec. This     *
 * structure contains IDASpgmr solver-specific data.              *
 *                                                                *
 ******************************************************************/


typedef struct {

  int  g_gstype;      /* type of Gram-Schmidt orthogonalization       */
  real g_sqrtN;       /* sqrt(N)                                      */
  int  g_maxl;        /* maxl = maximum dimension of the Krylov space */
  int  g_maxrs;       /* maxrs = max. number of GMRES restarts        */
  real g_eplifac;     /* eplifac = optional linear convergence factor */
  real g_dqincfac;    /* dqincfac = optional increment factor in Jv   */
  real g_epslin;      /* SpgrmSolve tolerance parameter               */

  int g_resflag;       /* flag from last res call                     */
  long int g_npe;      /* npe = total number of precond calls         */   
  long int g_nli;      /* nli = total number of linear iterations     */
  long int g_nps;      /* nps = total number of psolve calls          */
  long int g_ncfl;     /* ncfl = total number of convergence failures */

  long int g_nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int g_nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int g_nli0;     /* nli0 = saved nli (for performance monitor)   */   
  long int g_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int g_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int g_nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector g_ytemp;    /* temp vector used by IDAAtimesDQ            */ 
  N_Vector g_yptemp;   /* temp vector used by IDAAtimesDQ            */ 
  N_Vector g_xx;       /* temp vector used by IDASpgmrSolve          */
  N_Vector g_ycur;     /* IDA current y vector in Newton iteration   */
  N_Vector g_ypcur;    /* IDA current yp vector in Newton iteration  */
  N_Vector g_rcur;     /* rcur = F(tn, ycur, ypcur)                  */

  IDASpgmrPrecondFn g_precond; /* precond = user-supplied routine to   */
                               /* compute a preconditioner             */

  IDASpgmrPSolveFn g_psolve;   /* psolve = user-supplied routine to    */
			       /* solve preconditioner linear system   */

  void *g_pdata;               /* pdata passed to psolve and precond   */
  SpgmrMem g_spgmr_mem;        /* spgmr_mem is memory used by the      */
                               /* generic Spgmr solver                 */

} IDASpgmrMemRec, *IDASpgmrMem;


/* SensIDASPGMR linit, lsetup, lsolve, lperf, and lfree routines */

static int SensIDASpgmrInit(IDAMem ida_mem, boole *setupNonNull);

static int SensIDASpgmrSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
                          N_Vector resp, N_Vector tempv1,
                          N_Vector tempv2, N_Vector tempv3);

static int SensIDASpgmrSolve(IDAMem ida_mem, N_Vector bb, N_Vector ynow,
                          N_Vector ypnow, N_Vector rnow);

static int SensIDASpgmrPerf(IDAMem ida_mem, int perftask);

static int SensIDASpgmrFree(IDAMem ida_mem);


/* SensIDASPGMR Atimes and PSolve routines called by generic SPGMR solver */

static int SensIDASpgmrAtimesDQ(void *lin_mem, N_Vector v, N_Vector z);

static int SensIDASpgmrPSolve(void *lin_mem, N_Vector r, N_Vector z, int lr);


/* Readability Replacements */

#define uround  (ida_mem->ida_uround)
#define nst     (ida_mem->ida_nst)
#define tn      (ida_mem->ida_tn)
#define hh      (ida_mem->ida_hh)
#define cj      (ida_mem->ida_cj)
#define epsNewt (ida_mem->ida_epsNewt)
#define nre     (ida_mem->ida_nre)
#define res     (ida_mem->ida_res)
#define rdata   (ida_mem->ida_rdata)
#define ewt     (ida_mem->ida_ewt)
#define constraints (ida_mem->ida_constraints)
#define errfp   (ida_mem->ida_errfp)
#define iopt    (ida_mem->ida_iopt)
#define linit   (ida_mem->ida_linit)
#define lsetup  (ida_mem->ida_lsetup)
#define lsolve  (ida_mem->ida_lsolve)
#define lperf   (ida_mem->ida_lperf)
#define lfree   (ida_mem->ida_lfree)
#define lmem    (ida_mem->ida_lmem)
#define machenv (ida_mem->ida_machenv)
#define nni     (ida_mem->ida_nni)
#define ncfn    (ida_mem->ida_ncfn)

#define sqrtN   (IDASpgmr_mem->g_sqrtN)
#define epslin  (IDASpgmr_mem->g_epslin)
#define ytemp   (IDASpgmr_mem->g_ytemp)
#define yptemp  (IDASpgmr_mem->g_yptemp)
#define xx      (IDASpgmr_mem->g_xx)
#define ycur    (IDASpgmr_mem->g_ycur)
#define ypcur   (IDASpgmr_mem->g_ypcur)
#define rcur    (IDASpgmr_mem->g_rcur)
#define resflag (IDASpgmr_mem->g_resflag)
#define npe     (IDASpgmr_mem->g_npe)
#define nli     (IDASpgmr_mem->g_nli)
#define nps     (IDASpgmr_mem->g_nps)
#define ncfl    (IDASpgmr_mem->g_ncfl)
#define nst0    (IDASpgmr_mem->g_nst0)
#define nni0    (IDASpgmr_mem->g_nni0)
#define nli0    (IDASpgmr_mem->g_nli0)
#define ncfn0   (IDASpgmr_mem->g_ncfn0)
#define ncfl0   (IDASpgmr_mem->g_ncfl0)
#define nwarn   (IDASpgmr_mem->g_nwarn)

#define spgmr_mem (IDASpgmr_mem->g_spgmr_mem)


/*************** SensIDASpgmr *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the SensIDASPGMR linear solver module.  SensIDASpgmr sets
 the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and ida_lfree fields
 in (*IDA_mem) to be SensIDASpgmrInit, SensIDASpgmrSetup, SensIDASpgmrSolve,
 SensIDASpgmrPerf, and SensIDASpgmrFree, respectively.
 It allocates memory for a structure of type IDASpgmrMemRec and sets
 the ida_lmem field in (*IDA_mem) to the address of this structure.

 IDASpgmr sets the following fields in the IDASpgmrMemRec structure:

   g_gstype   = gstype
   g_maxl     = MIN(N,IDA_SPGMR_MAXL) if maxl <= 0,  else MIN(N,maxl)
   g_maxrs    = 0 if maxrs < 0,  MIN(5,N/g_maxl) if maxrs = 0, and
                MIN(maxrs,N/g_maxl) if maxrs > 0.
   g_eplifac  = 1.0 if eplifac = 0.0,  else eplifac
   g_dqincfac = 1.0 if dqincfac = 0.0,  else dqincfac
   g_pdata    = pdata
   g_precond  = precond
   g_psolve   = psolve

 SensIDASpgmr returns SUCCESS = 0 if successful, or 
          IDA_SPGMR_FAIL if IDA_mem is NULL or malloc failed, or
          IDA_SPGMR_BAD_ARG if the gstype argument is illegal.

**********************************************************************/

int SensIDASpgmr(void *IDA_mem, IDASpgmrPrecondFn precond, 
             IDASpgmrPSolveFn psolve, int gstype, int maxl, int maxrs,
             real eplifac, real dqincfac, void *pdata)

{
  IDAMem ida_mem;
  IDASpgmrMem IDASpgmr_mem;
  int maxl1, maxrs1;
  int N;
  SensData sdata;

  /* Return immediately if IDA_mem is NULL */
  ida_mem = (IDAMem) IDA_mem;
  if (ida_mem == NULL) return(IDA_SPGMR_FAIL);

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  /* Set five main function fields in ida_mem */
  linit  = SensIDASpgmrInit;
  lsetup = SensIDASpgmrSetup;
  lsolve = SensIDASpgmrSolve;
  lperf  = SensIDASpgmrPerf;
  lfree  = SensIDASpgmrFree;

  /* Get memory for IDASpgmrMemRec */
  lmem = IDASpgmr_mem = (IDASpgmrMem) malloc(sizeof(IDASpgmrMemRec));
  if (IDASpgmr_mem == NULL) return(IDA_SPGMR_FAIL);

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    fprintf(errfp, MSG_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
    return(IDA_SPGMR_BAD_ARG);
  }

  /* Set SPGMR parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPGMR_MAXL : maxl; maxl1 = MIN(maxl1, N);
  IDASpgmr_mem->g_gstype   = gstype;
  IDASpgmr_mem->g_maxl     = maxl1;
  maxrs1 = (maxrs == 0) ? IDA_SPGMR_MAXRS : maxrs;
  IDASpgmr_mem->g_maxrs    = (maxrs < 0) ? 0 : MIN(maxrs1, N/maxl1);
  IDASpgmr_mem->g_eplifac  = (eplifac == ZERO) ? ONE : eplifac;
  IDASpgmr_mem->g_dqincfac = (dqincfac == ZERO) ? ONE : dqincfac;
  IDASpgmr_mem->g_precond  = precond;
  IDASpgmr_mem->g_psolve   = psolve;
  IDASpgmr_mem->g_pdata    = pdata;


  /* Initialize counters and other optional outputs. */
  npe = nli = nps = ncfl = 0;
    
  if (iopt != NULL) {
    iopt[SPGMR_NPE] = npe;
    iopt[SPGMR_NLI] = nli;
    iopt[SPGMR_NPS] = nps;
    iopt[SPGMR_NCFL] = ncfl;
    iopt[SPGMR_LRW] = N*(maxl1 + 5) + maxl1*(maxl1 + 4) + 1;
    iopt[SPGMR_LIW] = 0;
  }

  return(SUCCESS);
}


/* Additional readability Replacements */

#define gstype   (IDASpgmr_mem->g_gstype)
#define maxl     (IDASpgmr_mem->g_maxl)
#define maxrs    (IDASpgmr_mem->g_maxrs)
#define eplifac  (IDASpgmr_mem->g_eplifac)
#define dqincfac (IDASpgmr_mem->g_dqincfac)
#define psolve   (IDASpgmr_mem->g_psolve)
#define precond  (IDASpgmr_mem->g_precond)
#define pdata    (IDASpgmr_mem->g_pdata)


/*************** SensIDASpgmrInit *****************************************

 This routine initializes remaining memory specific to the IDASPGMR 
 linear solver.  If any memory request fails, all memory previously
 allocated is freed, and an error message printed, before returning.

**********************************************************************/

static int SensIDASpgmrInit(IDAMem ida_mem, boole *setupNonNull)
{
  IDASpgmrMem IDASpgmr_mem;
  int N;
  SensData sdata;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  IDASpgmr_mem = (IDASpgmrMem) lmem;

  /* Print error message and return if IDASpgmr_mem is NULL */  
  if (IDASpgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }

  /* Allocate memory for ytemp, yptemp, and xx */
  ytemp = N_VNew(N, machenv);
  if (ytemp == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }
  yptemp = N_VNew(N, machenv);
  if (yptemp == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    return(LINIT_ERR);
  }
  xx = N_VNew(N, machenv);
  if (xx == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    N_VFree(yptemp);
    return(LINIT_ERR);
  }

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(N, maxl, machenv);
  if (spgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    N_VFree(yptemp);
    N_VFree(xx);
    return(LINIT_ERR);
  }

  /* Initialize sqrtN. */
  sqrtN = RSqrt(N);

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  *setupNonNull = (psolve != NULL) && (precond != NULL);

  return(LINIT_OK);
}

/*************** SensIDASpgmrSetup ****************************************

 This routine calls the user's preconditioner setup routine precond,
 and updates the counter npe.
 The return value is either
     SUCCESS = 0           if successful,
     LSETUP_ERROR_RECVR    if precond failed recoverably, or
     LSETUP_ERROR_NONRECVR if precond failed unrecoverably.

**********************************************************************/

static int SensIDASpgmrSetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp,
                         N_Vector resp, N_Vector tempv1, N_Vector tempv2,
                         N_Vector tempv3)
{
  int retval;
  IDASpgmrMem IDASpgmr_mem;
  int N;
  SensData sdata;
  N_Vector *yypsub, *yppsub, *respsub, *ewtsub;
  N_Vector *constraintssub, pconstraints;
  N_Vector *tempv1sub, *tempv2sub, *tempv3sub;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  IDASpgmr_mem = (IDASpgmrMem) lmem;

  /* Call user setup routine precond and update counter npe. */

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
   
  retval = precond(N, tn, yypsub[0], yppsub[0], respsub[0], cj, 
		   sdata->res_ptr, sdata->res_data, pdata, ewtsub[0], 
		   pconstraints, hh, uround, &nre, tempv1sub[0], 
		   tempv2sub[0], tempv3sub[0]);
  npe++;
  if (iopt != NULL) iopt[SPGMR_NPE] = npe;

  /* Return flag showing success or failure of precond. */
  if (retval < 0) return(LSETUP_ERROR_NONRECVR);
  if (retval > 0) return(LSETUP_ERROR_RECVR);
  return(SUCCESS);
}

/*************** SensIDASpgmrSolve ****************************************

 This routine handles the call to the generic SPGMR solver SpgmrSolve
 for the solution of the (1+Ns) linear systems Ax = b_i.

 The x-scaling and b-scaling arrays are both equal to ewt.

 We set the initial guess, x = 0, then call SpgmrSolve.  
 We copy the solution x into b_i, and update the counters nli, nps, ncfl.
 If SpgmrSolve returned nli_inc = 0 (hence x = 0), we take the SPGMR
 vtemp vector (= P_inverse F) as the correction vector instead.
 Finally, we set the return value according to the success of SpgmrSolve.

**********************************************************************/

static int SensIDASpgmrSolve(IDAMem ida_mem, N_Vector bb, N_Vector ynow,
			N_Vector ypnow, N_Vector rnow)
{
  IDASpgmrMem IDASpgmr_mem;
  int pretype, nli_inc, nps_inc, retval;
  real res_norm;
  int N;
  int i, Ns;
  N_Vector *ynowsub, *ypnowsub, *rnowsub;
  N_Vector *bbsub, *ewtsub;
  SensData sdata;

  sdata = (SensData) ida_mem->ida_rdata;
  ynowsub = N_VSUB(ynow);
  ypnowsub = N_VSUB(ypnow);
  rnowsub = N_VSUB(rnow);
  bbsub = N_VSUB(bb);
  ewtsub = N_VSUB(ewt);
  N = sdata->Ny;
  Ns = sdata->Ns;

  IDASpgmr_mem = (IDASpgmrMem) lmem;

  /* Set SpgmrSolve convergence test constant epslin, in terms of the
    Newton convergence test constant epsNewt and safety factors.  The factor 
    sqrt(N) assures that the GMRES convergence test is applied to the
    WRMS norm of the residual vector, rather than the weighted L2 norm. */
  epslin = sqrtN*eplifac*PT05*epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  ycur = ynowsub[0];
  ypcur = ypnowsub[0];
  rcur = rnowsub[0];

  /* Set SpgmrSolve inputs pretype and initial guess xx = 0. */  
  pretype = (psolve == NULL) ? NONE : LEFT;

  for (i = 0; i <= Ns; i++) {
    N_VConst(ZERO, xx);
  
    /* Call SpgmrSolve and copy xx to bb. */
    retval = SpgmrSolve(spgmr_mem, ida_mem, xx, bbsub[i], pretype, gstype, 
			epslin, maxrs, ida_mem, ewtsub[i], ewtsub[i], 
			SensIDASpgmrAtimesDQ, SensIDASpgmrPSolve, &res_norm, 
			&nli_inc, &nps_inc);

    if (nli_inc == 0) N_VScale(ONE, SPGMR_VTEMP(spgmr_mem), bbsub[i]);
    else N_VScale(ONE, xx, bbsub[i]);
    
    /* Increment counters nli, nps, and return if successful. */
    nli += nli_inc;
    nps += nps_inc;
    if (iopt != NULL) {
      iopt[SPGMR_NLI] = nli;
      iopt[SPGMR_NPS] = nps;
    }
    if (retval != SUCCESS) break;
  }

  if (retval == SUCCESS) return(SUCCESS);

  /* If not successful, increment ncfl and return appropriate flag. */
  ncfl++;
  if (iopt != NULL) iopt[SPGMR_NCFL] = ncfl;

  if (retval > 0)   return(LSOLVE_ERROR_RECVR);
  if (retval != -2) return(LSOLVE_ERROR_NONRECVR);
  if (resflag > 0)  return(LSOLVE_ERROR_RECVR);
                    return(LSOLVE_ERROR_NONRECVR);

}

/*************** SensIDASpgmrPerf *****************************************

 This routine handles performance monitoring specific to the IDASPGMR
 linear solver.  When perftask = 0, it saves values of various counters.
 When perftask = 1, it examines difference quotients in these counters
 and, depending on their values, it prints up to three warning messages.
 Messages are printed up to a maximum of 10 times.

**********************************************************************/

static int SensIDASpgmrPerf(IDAMem ida_mem, int perftask)
{
  IDASpgmrMem IDASpgmr_mem;
  real avdim, rcfn, rcfl;
  long int nstd, nnid;
  boole lavd, lcfn, lcfl;
  SensData sdata;
  int Ns;

  IDASpgmr_mem = (IDASpgmrMem) lmem;

  sdata = (SensData) ida_mem->ida_rdata;
  Ns = sdata->Ns;

  if (perftask == 0) {
    nst0 = nst;  nni0 = nni;  nli0 = nli;
    ncfn0 = ncfn;  ncfl0 = ncfl;  
    nwarn = 0;
    return;
  }

  nstd = nst - nst0;  nnid = nni - nni0;
  if (nstd == 0 || nnid == 0) return;
  avdim = (nli - nli0)/( (real) nnid * (1+Ns));
  rcfn = (ncfn - ncfn0)/( (real) nstd);
  rcfl = (ncfl - ncfl0)/( (real) nnid * (1+Ns));
  lavd = (avdim > ( (real) maxl ) );
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lavd || lcfn || lcfl)) return;
  nwarn++;
  if (nwarn > 10) return;
  if (lavd) fprintf(errfp, MSG_AVD_WARN, tn, avdim);
  if (lcfn) fprintf(errfp, MSG_CFN_WARN, tn, rcfn);
  if (lcfl) fprintf(errfp, MSG_CFL_WARN, tn, rcfl);
}

/*************** SensIDASpgmrFree *****************************************

 This routine frees memory specific to the IDASPGMR linear solver.

**********************************************************************/

static int SensIDASpgmrFree(IDAMem ida_mem)
{
  IDASpgmrMem IDASpgmr_mem;

  IDASpgmr_mem = (IDASpgmrMem) lmem;
  
  N_VFree(ytemp);
  N_VFree(yptemp);
  N_VFree(xx);
  SpgmrFree(spgmr_mem);
  free(lmem);

}

/*************** SensIDASpgmrAtimesDQ *************************************

 This routine generates the matrix-vector product z = Jv, where
 J is the system Jacobian, by using a difference quotient approximation.
 The approximation is 
      Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  where
        y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
        sigma = sqrt(N)*dqincfac.
 The return value from the call to res is saved in order to set the
 return flag from SensIDASpgmrSolve.

**********************************************************************/

static int SensIDASpgmrAtimesDQ(void *idamem, N_Vector v, N_Vector z)
{
  IDAMem   ida_mem;
  IDASpgmrMem IDASpgmr_mem;
  real sig, siginv;
  int ires;
  int N;
  SensData sdata;

  ida_mem = (IDAMem) idamem;
  IDASpgmr_mem = (IDASpgmrMem) lmem;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;

  sig = sqrtN*dqincfac;

  /* Set ytemp = ycur + sig*v, yptemp = ypcur + cj*sig*v. */
  N_VLinearSum(sig, v, ONE, ycur, ytemp);
  N_VLinearSum(cj*sig, v, ONE, ypcur, yptemp);

  /* Call res for z = F(t, ytemp, yptemp), and return if it failed. */
  ires = sdata->res_ptr(N, tn, ytemp, yptemp, z, sdata->res_data); 
  nre++;
  resflag = ires;
  if (ires != 0) return(ires);

  /* Set z to [z - rcur]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, z, -siginv, rcur, z);

  return(SUCCESS);

}

/*************** SensIDASpgmrPSolve ***************************************

 This routine interfaces between the generic SpgmrSolve routine and
 the user's psolve routine.  It passes to psolve all required state 
 information from ida_mem.  Its return value is the same as that
 returned by psolve.  Note that the generic SPGMR solver guarantees
 that SensIDASpgmrPSolve will not be called in the case psolve = NULL.

**********************************************************************/

static int SensIDASpgmrPSolve(void *idamem, N_Vector r, N_Vector z, int lr)
{
  IDAMem   ida_mem;
  IDASpgmrMem IDASpgmr_mem;
  int retval;
  int N;
  N_Vector *ewtsub;
  SensData sdata;

  ida_mem = (IDAMem) idamem;
  IDASpgmr_mem = (IDASpgmrMem) lmem;

  sdata = (SensData) ida_mem->ida_rdata;
  N = sdata->Ny;
  ewtsub = N_VSUB(ewt);

  retval = psolve(N, tn, ycur, ypcur, rcur, cj, sdata->res_ptr, 
		  sdata->res_data, pdata, ewtsub[0], epslin, r, z, 
		  &nre, ytemp);
  /* This call is counted in nps within the SensIDASpgmrSolve routine */

  return(retval);

}
