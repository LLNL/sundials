/*******************************************************************
 *                                                                 *
 * File          : cvspgmr.c                                       *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh and           * 
 *                 Radu Serban @ LLNL                              *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the implementation file for the CVODE scaled,           *
 * preconditioned GMRES linear solver, CVSPGMR.                    *
 *                                                                 *
 *******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "cvspgmr.h"
#include "cvode.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "iterativ.h"
#include "spgmr.h"


/* Error Messages */

#define CVSPGMR            "CVSpgmr/CVReInitSpgmr-- "

#define MSG_CVMEM_NULL     CVSPGMR "CVode Memory is NULL.\n\n"

#define MSG_MEM_FAIL       CVSPGMR "A memory request failed.\n\n"

#define MSG_BAD_PRETYPE_1  CVSPGMR "pretype=%d illegal.\n"
#define MSG_BAD_PRETYPE_2  "The legal values are NONE=%d, LEFT=%d, "
#define MSG_BAD_PRETYPE_3  "RIGHT=%d, and BOTH=%d.\n\n"
#define MSG_BAD_PRETYPE    MSG_BAD_PRETYPE_1 MSG_BAD_PRETYPE_2 MSG_BAD_PRETYPE_3

#define MSG_PSOLVE_REQ_1   CVSPGMR "pretype!=NONE, but PSOLVE=NULL is "
#define MSG_PSOLVE_REQ_2   "illegal.\n\n"
#define MSG_PSOLVE_REQ     MSG_PSOLVE_REQ_1 MSG_PSOLVE_REQ_2

#define MSG_BAD_GSTYPE_1   CVSPGMR "gstype=%d illegal.\n"
#define MSG_BAD_GSTYPE_2   "The legal values are MODIFIED_GS=%d and "
#define MSG_BAD_GSTYPE_3   "CLASSICAL_GS=%d.\n\n"
#define MSG_BAD_GSTYPE     MSG_BAD_GSTYPE_1 MSG_BAD_GSTYPE_2 MSG_BAD_GSTYPE_3

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/******************************************************************
 *                                                                *           
 * Types : CVSpgmrMemRec, CVSpgmrMem                              *
 *----------------------------------------------------------------*
 * The type CVSpgmrMem is pointer to a CVSpgmrMemRec. This        *
 * structure contains CVSpgmr solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {

  int  g_pretype;     /* type of preconditioning                      */
  int  g_gstype;      /* type of Gram-Schmidt orthogonalization       */
  realtype g_sqrtN;   /* sqrt(N)                                      */
  realtype g_delt;    /* delt = user specified or DELT_DEFAULT        */
  realtype g_deltar;  /* deltar = delt * tq4                          */
  realtype g_delta;   /* delta = deltar * sqrtN                       */
  int  g_maxl;        /* maxl = maximum dimension of the Krylov space */

  long int g_nstlpre;  /* value of nst at the last precond call       */     
  long int g_npe;      /* npe = total number of precond calls         */   
  long int g_nli;      /* nli = total number of linear iterations     */
  long int g_nps;      /* nps = total number of psolve calls          */
  long int g_ncfl;     /* ncfl = total number of convergence failures */

  N_Vector g_ytemp;    /* temp vector passed to jtimes and psolve     */ 
  N_Vector g_x;        /* temp vector used by CVSpgmrSolve            */
  N_Vector g_ycur;     /* CVODE current y vector in Newton Iteration  */
  N_Vector g_fcur;     /* fcur = f(tn, ycur)                          */

  CVSpgmrPrecondFn g_precond; /* precond = user-supplied routine to   */
                              /* compute a preconditioner             */

  CVSpgmrPSolveFn g_psolve;   /* psolve = user-supplied routine to    */
                              /* solve preconditioner linear system   */

  void *g_P_data;            /* P_data passed to psolve and precond   */
  SpgmrMem g_spgmr_mem;      /* spgmr_mem is memory used by the       */
                             /* generic Spgmr solver                  */

  CVSpgmrJtimesFn g_jtimes;  /* jtimes = Jacobian * vector routine    */
                             /*          to be called                 */
  void *g_j_data;            /* j_data is passed to jtimes            */

} CVSpgmrMemRec, *CVSpgmrMem;


/* CVSPGMR linit, lsetup, lsolve, and lfree routines */

static int CVSpgmrInit(CVodeMem cv_mem);

static int CVSpgmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3);

static int CVSpgmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector ycur,
                        N_Vector fcur);

static void CVSpgmrFree(CVodeMem cv_mem);

/* CVSPGMR Atimes and PSolve routines called by generic SPGMR solver */

static int CVSpgmrAtimes(void *cv_mem, N_Vector v, N_Vector z);

static int CVSpgmrPSolve(void *cv_mem, N_Vector r, N_Vector z, int lr);

/* CVSPGMR difference quotient routine for J*v */

static int CVSpgmrDQJtimes(integertype N, N_Vector v, N_Vector Jv, RhsFn f,
                       void *f_data, realtype tn, N_Vector y, N_Vector fy,
                       realtype vnrm, N_Vector ewt, realtype h, realtype uround,
                       void *jac_data, long int *nfePtr, N_Vector ytemp);


/*************** CVSpgmrDQJtimes *************************************

 This routine generates a difference quotient approximation to
 the Jacobian times vector f_y(t,y) * v. The approximation is 
 Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 input, i.e. the WRMS norm of v/vnrm is 1.

**********************************************************************/

static int CVSpgmrDQJtimes(integertype N, N_Vector v, N_Vector Jv, RhsFn f, 
                       void *f_data, realtype tn, N_Vector y, N_Vector fy,
                       realtype vnrm, N_Vector ewt, realtype h, realtype uround,
                       void *jac_data, long int *nfePtr, N_Vector work)
{
  N_Vector ytemp;

  /* Rename work for readability */
  ytemp = work;

  /* Set ytemp = y + (1/vnrm) v */
  N_VLinearSum(ONE/vnrm, v, ONE, y, ytemp);

  /* Set Jv = f(tn, ytemp) */
  f(N, tn, ytemp, Jv, f_data); 
  (*nfePtr)++;
  
  /* Replace Jv by vnrm*(Jv - fy) */
  N_VLinearSum(ONE, Jv, -ONE, fy, Jv);
  N_VScale(vnrm, Jv, Jv);

  return(0);
}

/**********************************************************************/

/* Readability Replacements */

#define N       (cv_mem->cv_N)      
#define uround  (cv_mem->cv_uround)
#define tq      (cv_mem->cv_tq)
#define nst     (cv_mem->cv_nst)
#define tn      (cv_mem->cv_tn)
#define h       (cv_mem->cv_h)
#define gamma   (cv_mem->cv_gamma)
#define gammap  (cv_mem->cv_gammap)   
#define nfe     (cv_mem->cv_nfe)
#define f       (cv_mem->cv_f)
#define f_data  (cv_mem->cv_f_data)
#define ewt     (cv_mem->cv_ewt)
#define errfp   (cv_mem->cv_errfp)
#define mnewt   (cv_mem->cv_mnewt)
#define iopt    (cv_mem->cv_iopt)
#define ropt    (cv_mem->cv_ropt)
#define linit   (cv_mem->cv_linit)
#define lsetup  (cv_mem->cv_lsetup)
#define lsolve  (cv_mem->cv_lsolve)
#define lfree   (cv_mem->cv_lfree)
#define lmem    (cv_mem->cv_lmem)
#define machenv (cv_mem->cv_machenv)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define sqrtN   (cvspgmr_mem->g_sqrtN)   
#define ytemp   (cvspgmr_mem->g_ytemp)
#define x       (cvspgmr_mem->g_x)
#define ycur    (cvspgmr_mem->g_ycur)
#define fcur    (cvspgmr_mem->g_fcur)
#define delta   (cvspgmr_mem->g_delta)
#define deltar  (cvspgmr_mem->g_deltar)
#define npe     (cvspgmr_mem->g_npe)
#define nli     (cvspgmr_mem->g_nli)
#define nps     (cvspgmr_mem->g_nps)
#define ncfl    (cvspgmr_mem->g_ncfl)
#define nstlpre (cvspgmr_mem->g_nstlpre)
#define spgmr_mem (cvspgmr_mem->g_spgmr_mem)


/*************** CVSpgmr *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the Spgmr linear solver module. CVSpgmr first
 calls the existing lfree routine if this is not NULL.  It then sets
 the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 to be CVSpgmrInit, CVSpgmrSetup, CVSpgmrSolve, and CVSpgmrFree,
 respectively.  It allocates memory for a structure of type
 CVSpgmrMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure.  It sets setupNonNull in (*cvode_mem),
 and sets the following fields in the CVSpgmrMemRec structure:
   g_pretype = pretype                                       
   g_gstype  = gstype                                       
   g_maxl    = MIN(N,CVSPGMR_MAXL)  if maxl <= 0             
             = maxl                 if maxl > 0              
   g_delt    = CVSPGMR_DELT if delt == 0.0                     
             = delt         if delt != 0.0                     
   g_P_data  = P_data                                        
   g_precond = precond                                       
   g_psolve  = psolve                                        
   g_jtimes  = input parameter jtimes  if jtimes != NULL
             = CVSpgmrDQJtimes         otherwise
   g_j_data  = input parameter jac_data
 Finally, CVSpgmr allocates memory for ytemp and x, and calls
 SpgmrMalloc to allocate memory for the Spgmr solver.  The CVSpgmr
 return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

**********************************************************************/

int CVSpgmr(void *cvode_mem, int pretype, int gstype, int maxl, realtype delt,
             CVSpgmrPrecondFn precond, CVSpgmrPSolveFn psolve, void *P_data,
             CVSpgmrJtimesFn jtimes, void *jac_data)

{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem == NULL) {                /* CVode reports this error */
    fprintf(errfp, MSG_CVMEM_NULL);
    return(LMEM_FAIL);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = CVSpgmrInit;
  lsetup = CVSpgmrSetup;
  lsolve = CVSpgmrSolve;
  lfree  = CVSpgmrFree;

  /* Get memory for CVSpgmrMemRec */
  lmem = cvspgmr_mem = (CVSpgmrMem) malloc(sizeof(CVSpgmrMemRec));
  if (cvspgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }

  /* Set Spgmr parameters that have been passed in call sequence */
  cvspgmr_mem->g_pretype    = pretype;
  cvspgmr_mem->g_gstype     = gstype;
  mxl = cvspgmr_mem->g_maxl = (maxl <= 0) ? MIN(CVSPGMR_MAXL, N) : maxl;
  cvspgmr_mem->g_delt       = (delt == ZERO) ? CVSPGMR_DELT : delt;
  cvspgmr_mem->g_P_data     = P_data;
  cvspgmr_mem->g_precond    = precond;
  cvspgmr_mem->g_psolve     = psolve;

  /* Set Jacobian times vector routine to user's jtimes or CVSpgmrDQJtimes */
  if(jtimes == NULL) {
    cvspgmr_mem->g_jtimes = CVSpgmrDQJtimes;
  } else {
    cvspgmr_mem->g_jtimes = jtimes;
  }

  /* Set Jacobian data */
  cvspgmr_mem->g_j_data = jac_data;

  /* Check for legal pretype, precond, and psolve */ 
  if ((pretype != NONE) && (pretype != LEFT) &&
      (pretype != RIGHT) && (pretype != BOTH)) {
    fprintf(errfp, MSG_BAD_PRETYPE, pretype, NONE, LEFT, RIGHT, BOTH);
    return(LIN_ILL_INPUT);
  }
  if ((pretype != NONE) && (psolve == NULL)) {
    fprintf(errfp, MSG_PSOLVE_REQ);
    return(LIN_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    fprintf(errfp, MSG_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
    return(LIN_ILL_INPUT);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning (pretype != NONE)
     and there is a preconditioning setup phase (precond != NULL)          */
  setupNonNull = (pretype != NONE) && (precond != NULL);

  /* Allocate memory for ytemp and x */
  ytemp = N_VNew(N, machenv);
  if (ytemp == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }
  x = N_VNew(N, machenv);
  if (x == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    return(LMEM_FAIL);
  }

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(N, mxl, machenv);
  if (spgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    N_VFree(x);
    return(LMEM_FAIL);
  }

  return(SUCCESS);
}


/*************** CVReInitSpgmr****************************************

 This routine resets the link between the main CVODE module and the
 Spgmr linear solver module CVSPGMR.  No memory freeing or allocation
 operations are done, as the existing linear solver memory is assumed
 sufficient.  All other initializations are the same as in CVSpgmr.
 The return value is SUCCESS=0, LMEM_FAIL=-1, or LIN_ILL_INPUT=-2.

**********************************************************************/

int CVReInitSpgmr(void *cvode_mem, int pretype, int gstype, int maxl,
            realtype delt, CVSpgmrPrecondFn precond, CVSpgmrPSolveFn psolve,
            void *P_data, CVSpgmrJtimesFn jtimes, void *jac_data)

{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem == NULL) {                /* CVode reports this error */
    fprintf(errfp, MSG_CVMEM_NULL);
    return(LMEM_FAIL);
  }

  /* Set four main function fields in cv_mem */
  linit  = CVSpgmrInit;
  lsetup = CVSpgmrSetup;
  lsolve = CVSpgmrSolve;
  lfree  = CVSpgmrFree;

  cvspgmr_mem = lmem;   /* Use existing linear solver memory pointer */

  /* Set Spgmr parameters that have been passed in call sequence */
  cvspgmr_mem->g_pretype    = pretype;
  cvspgmr_mem->g_gstype     = gstype;
  mxl = cvspgmr_mem->g_maxl = (maxl <= 0) ? MIN(CVSPGMR_MAXL, N) : maxl;
  cvspgmr_mem->g_delt       = (delt == ZERO) ? CVSPGMR_DELT : delt;
  cvspgmr_mem->g_P_data     = P_data;
  cvspgmr_mem->g_precond    = precond;
  cvspgmr_mem->g_psolve     = psolve;

  /* Set Jacobian times vector routine to user's jtimes or CVSpgmrDQJtimes */
  if(jtimes == NULL) {
    cvspgmr_mem->g_jtimes = CVSpgmrDQJtimes;
  } else {
    cvspgmr_mem->g_jtimes = jtimes;
  }

  /* Set Jacobian data */
  cvspgmr_mem->g_j_data = jac_data;

  /* Check for legal pretype, precond, and psolve */ 
  if ((pretype != NONE) && (pretype != LEFT) &&
      (pretype != RIGHT) && (pretype != BOTH)) {
    fprintf(errfp, MSG_BAD_PRETYPE, pretype, NONE, LEFT, RIGHT, BOTH);
    return(LIN_ILL_INPUT);
  }
  if ((pretype != NONE) && (psolve == NULL)) {
    fprintf(errfp, MSG_PSOLVE_REQ);
    return(LIN_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    fprintf(errfp, MSG_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
    return(LIN_ILL_INPUT);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning (pretype != NONE)
     and there is a preconditioning setup phase (precond != NULL)          */
  setupNonNull = (pretype != NONE) && (precond != NULL);

  return(SUCCESS);
}


/* Additional readability Replacements */

#define pretype (cvspgmr_mem->g_pretype)
#define gstype  (cvspgmr_mem->g_gstype)
#define delt    (cvspgmr_mem->g_delt)
#define maxl    (cvspgmr_mem->g_maxl)
#define psolve  (cvspgmr_mem->g_psolve)
#define precond (cvspgmr_mem->g_precond)
#define P_data  (cvspgmr_mem->g_P_data)
#define jtimes  (cvspgmr_mem->g_jtimes)
#define j_data  (cvspgmr_mem->g_j_data)


/*************** CVSpgmrInit *****************************************

 This routine does remaining initializations specific to the Spgmr 
 linear solver.

**********************************************************************/

static int CVSpgmrInit(CVodeMem cv_mem)
{
  CVSpgmrMem cvspgmr_mem;
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Initialize sqrtN and counters, and set workspace lengths */

  sqrtN = RSqrt(N);
  npe = nli = nps = ncfl = nstlpre = 0;
    
  if (iopt != NULL) {
    iopt[SPGMR_NPE] = npe;
    iopt[SPGMR_NLI] = nli;
    iopt[SPGMR_NPS] = nps;
    iopt[SPGMR_NCFL] = ncfl;
    iopt[SPGMR_LRW] = N*(maxl + 5) + maxl*(maxl + 4) + 1;
    iopt[SPGMR_LIW] = 0;
  }

  return(LINIT_OK);
}

/*************** CVSpgmrSetup ****************************************

 This routine does the setup operations for the Spgmr linear solver.
 It makes a decision as to whether or not to signal for re-evaluation
 of Jacobian data in the precond routine, based on various state
 variables, then it calls precond.  If we signal for re-evaluation,
 then we reset jcur = *jcurPtr to TRUE, regardless of the precond output.
 In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.

**********************************************************************/

static int CVSpgmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  ier;
  CVSpgmrMem cvspgmr_mem;

  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlpre + CVSPGMR_MSBPRE) ||
    ((convfail == FAIL_BAD_J) && (dgamma < CVSPGMR_DGMAX)) ||
    (convfail == FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call precond routine and possibly reset jcur */
  ier = precond(N, tn, ypred, fpred, jok, jcurPtr, gamma, ewt, h,
                uround, &nfe, P_data, vtemp1, vtemp2, vtemp3);
  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    npe++;
    nstlpre = nst;
  }

  /* Set npe, and return the same value ier that precond returned */
  if (iopt != NULL) iopt[SPGMR_NPE] = npe;
  return(ier);
}

/*************** CVSpgmrSolve ****************************************

 This routine handles the call to the generic solver SpgmrSolve
 for the solution of the linear system Ax = b with the SPGMR method,
 without restarts.  The solution x is returned in the vector b.

 If the WRMS norm of b is small, we return x = b (if this is the first
 Newton iteration) or x = 0 (if a later Newton iteration).

 Otherwise, we set the tolerance parameter and initial guess (x = 0),
 call SpgmrSolve, and copy the solution x into b.  The x-scaling and
 b-scaling arrays are both equal to ewt, and no restarts are allowed.

 The counters nli, nps, and ncfl are incremented, and the return value
 is set according to the success of SpgmrSolve.  The success flag is
 returned if SpgmrSolve converged, or if this is the first Newton
 iteration and the residual norm was reduced below its initial value.

**********************************************************************/

static int CVSpgmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector ynow,
                        N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSpgmrMem cvspgmr_mem;
  int nli_inc, nps_inc, ier;
  
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = delt*tq[4]; 

  bnorm = N_VWrmsNorm(b, ewt);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur = ynow;
  fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SpgmrSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SpgmrSolve and copy x to b */
  ier = SpgmrSolve(spgmr_mem, cv_mem, x, b, pretype, gstype, delta, 0,
                   cv_mem, ewt, ewt, CVSpgmrAtimes, CVSpgmrPSolve,
                   &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (iopt != NULL) {
    iopt[SPGMR_NLI] = nli;
    iopt[SPGMR_NPS] = nps;
  }  
  if (ier != 0) { 
    ncfl++;
    if (iopt != NULL) iopt[SPGMR_NCFL] = ncfl;
  }

  /* Set return value to -1, 0, or 1 */
  if (ier < 0) return(-1);  
  if ((ier == SPGMR_SUCCESS) || 
      ((ier == SPGMR_RES_REDUCED) && (mnewt == 0)))
    return(0);
  return(1);  
}

/*************** CVSpgmrFree *****************************************

 This routine frees memory specific to the Spgmr linear solver.

**********************************************************************/

static void CVSpgmrFree(CVodeMem cv_mem)
{
  CVSpgmrMem cvspgmr_mem;

  cvspgmr_mem = (CVSpgmrMem) lmem;
  
  N_VFree(ytemp);
  N_VFree(x);
  SpgmrFree(spgmr_mem);
  free(cvspgmr_mem);
}

/***************** CVSpgmrAtimes *************************************

 This routine generates the matrix-vector product z = Mv, where
 M = I - gamma*J. The product J*v is obtained by calling the jtimes 
 routine. It is then scaled by -gamma and added to v to obtain M*v.
 The return value is the same as the value returned by jtimes --
 0 if successful, nonzero otherwise.

**********************************************************************/

static int CVSpgmrAtimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  realtype rho;
  CVodeMem   cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int jtflag;

  cv_mem = (CVodeMem) cvode_mem;
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* If rho = norm(v) is 0, return z = 0 */
  rho = N_VWrmsNorm(v, ewt);
  if (rho == ZERO) {
    N_VConst(ZERO, z);
    return(0);
  }
  
  jtflag = jtimes(N, v, z, f, f_data, tn, ycur, fcur, rho, ewt, h,
                  uround, j_data, &nfe, ytemp);
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*************** CVSpgmrPSolve ***************************************

 This routine interfaces between the generic SpgmrSolve routine and
 the user's psolve routine.  It passes to psolve all required state 
 information from cvode_mem.  Its return value is the same as that
 returned by psolve. Note that the generic SPGMR solver guarantees
 that CVSpgmrPSolve will not be called in the case in which
 preconditioning is not done. This is the only case in which the
 user's psolve routine is allowed to be NULL.

**********************************************************************/

static int CVSpgmrPSolve(void *cvode_mem, N_Vector r, N_Vector z, int lr)
{
  CVodeMem   cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int ier;

  cv_mem = (CVodeMem) cvode_mem;
  cvspgmr_mem = (CVSpgmrMem)lmem;

  ier = psolve(N, tn, ycur, fcur, ytemp, gamma, ewt, delta, &nfe, r,
               lr, P_data, z);
  /* This call is counted in nps within the CVSpgmrSolve routine */

  return(ier);     
}
