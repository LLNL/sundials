/*******************************************************************
 *                                                                 *
 * File          : kinspgmr.c                                      *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 4 August 2003                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/kinsol/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for the KINSOL scaled,          *
 * preconditioned GMRES linear solver, KINSPGMR.                   *
 *                                                                 *
 *******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "kinspgmr.h"
#include "kinsol.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "iterative.h"
#include "spgmr.h"


/* Error Messages */

#define KINSPGMR      "KINSpgmr-- "

#define MSG_KINMEM_NULL KINSPGMR "KINSOL memory is NULL.\n\n"

#define MSG_MEM_FAIL    KINSPGMR "A memory request failed.\n\n"

#define MSG_SETGET_KINMEM_NULL "KINSpgmrSet*/KINSpgmrGet*-- KIN memory is NULL. \n\n"

#define MSG_SETGET_LMEM_NULL   "KINSpgmrSet*/KINSpgmrGet*-- KINSPGMR memory is NULL. \n\n"

#define MSG_KINS_NEG_MAXRS     "KINSpgmrSetMaxRestarts-- maxrs<0 illegal. \n\n"

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* KINSpgmr linit, lsetup, lsolve, and lfree routines */

static int KINSpgmrInit(KINMem kin_mem);

static int KINSpgmrSetup(KINMem kin_mem);

static int KINSpgmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb,
                         realtype *res_norm);

static int KINSpgmrFree(KINMem kin_mem);

/* KINSpgmr Atimes and PSolve routines called by generic SPGMR solver */

static int KINSpgmrAtimes(void *kinsol_mem, N_Vector v, N_Vector z);

static int KINSpgmrPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

static int KINSpgmrDQJtimes(N_Vector v, N_Vector Jv,
                            N_Vector u, booleantype *new_u, 
                            void *jac_data);


/* Readability Replacements */

#define lrw1         (kin_mem->kin_lrw1)
#define liw1         (kin_mem->kin_liw1)
#define uround       (kin_mem->kin_uround)
#define nfe          (kin_mem->kin_nfe)
#define nni          (kin_mem->kin_nni)
#define nnilpre      (kin_mem->kin_nnilpre)
#define func         (kin_mem->kin_func)
#define f_data       (kin_mem->kin_f_data)
#define printfl      (kin_mem->kin_printfl)
#define linit        (kin_mem->kin_linit)
#define lsetup       (kin_mem->kin_lsetup)
#define lsolve       (kin_mem->kin_lsolve)
#define lfree        (kin_mem->kin_lfree)
#define lmem         (kin_mem->kin_lmem)
#define uu           (kin_mem->kin_uu)
#define fval         (kin_mem->kin_fval)
#define uscale       (kin_mem->kin_uscale)
#define fscale       (kin_mem->kin_fscale)
#define sqrt_relfunc (kin_mem->kin_sqrt_relfunc)
#define precondflag  (kin_mem->kin_precondflag)
#define precondcurrent (kin_mem->kin_precondcurrent)
#define eps          (kin_mem->kin_eps)
#define sJpnorm      (kin_mem->kin_sJpnorm)
#define sfdotJp      (kin_mem->kin_sfdotJp)
#define errfp        (kin_mem->kin_errfp)
#define setupNonNull (kin_mem->kin_setupNonNull)
#define vtemp1       (kin_mem->kin_vtemp1)
#define vtemp2       (kin_mem->kin_vtemp2)
#define nvspec       (kin_mem->kin_nvspec)

#define pretype (kinspgmr_mem->g_pretype)
#define gstype  (kinspgmr_mem->g_gstype)
#define nli     (kinspgmr_mem->g_nli)
#define npe     (kinspgmr_mem->g_npe)
#define nps     (kinspgmr_mem->g_nps)
#define ncfl    (kinspgmr_mem->g_ncfl)
#define njtimes (kinspgmr_mem->g_njtimes)
#define nfeSG   (kinspgmr_mem->g_nfeSG)
#define new_uu  (kinspgmr_mem->g_new_uu)

#define spgmr_mem (kinspgmr_mem->g_spgmr_mem)

/*************** KINSpgmr *********************************************

 This routine allocates and initializes the memory record and sets various 
 function fields specific to the Spgmr linear solver module. KINSpgmr sets 
 the kin_linit, kin_lsetup, kin_lsolve, and kin_lfree fields in (*kinsol_mem)
 to be KINSpgmrInit, KINSpgmrSetup, KINSpgmrSolve, and KINSpgmrFree,
 respectively. It sets fields in the kinspgmr_mem block for the preconditioner
 setup and solve routines based on input parameters precondset and psolve. It
 also sets the kinspgmr_mem entry element for the optional user-supplied
 Atimes (Jacobian J * v) routine userAtimes. It allocates memory for a 
 structure of type KINSpgmrMemRec and sets the kin_lmem field in 
 (*kinsol_mem) to the address of this structure. It also calls SpgmrMalloc
 to allocate memory for the module SPGMR. In summary, KINSpgmr sets 
 the following fields in the KINSpgmrMemRec structure: 
                          
   pretype   = RIGHT, if the PrecondSolve routine is provided, else NONE
   gstype    = MODIFIED_GS
   g_maxl    = MIN(Neq,KINSPGMR_MAXL)  if maxl <= 0             
             = maxl                    if maxl > 0   
   g_maxlrst = maxlrst           
   g_precond = precondset
   g_psolve  = psolve    
   g_userAtimes = userAtimes                                    
   g_P_data  = P_data                                        

**********************************************************************/

int KINSpgmr(void *kinmem, int maxl)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;
  int maxl1;

  if (kinmem == NULL){
     fprintf(stdout, MSG_KINMEM_NULL);
     return(KIN_MEM_NULL);  
  }
  kin_mem = (KINMem) kinmem;

  /* Set four main function fields in kin_mem */
  linit  = KINSpgmrInit; 
  lsetup = KINSpgmrSetup;
  lsolve = KINSpgmrSolve;
  lfree  = KINSpgmrFree;

  /* Get memory for KINSpgmrMemRec */
  kinspgmr_mem = (KINSpgmrMem) malloc(sizeof(KINSpgmrMemRec));
  if (kinspgmr_mem == NULL){
     fprintf(errfp, MSG_MEM_FAIL);
     return(KINSPGMR_MEM_FAIL);  
  }

  /* Set SPGMR parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? KINSPGMR_MAXL : maxl;
  kinspgmr_mem->g_maxl = maxl1;  

  /* Set default values for the rest of the Spgmr parameters */
  kinspgmr_mem->g_pretype  = NONE;
  kinspgmr_mem->g_gstype   = MODIFIED_GS;
  kinspgmr_mem->g_maxlrst  = 0;
  kinspgmr_mem->g_P_data   = NULL;
  kinspgmr_mem->g_pset     = NULL;
  kinspgmr_mem->g_psolve   = NULL;
  kinspgmr_mem->g_jtimes   = NULL;

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(maxl1, nvspec);
  if (spgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    lmem = NULL; /* Set lmem to NULL and free that memory, as a flag to a
                    later inadvertent KINSol call, that SpgmrMalloc failed. */
    free(lmem);
    return(SPGMR_MEM_FAIL);
  }

  /* Attach linear solver memory to KINSOL memory */
  lmem = kinspgmr_mem;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrSetMaxRestarts(void *kinmem, int maxrs)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Check for legal maxrs */
  if (maxrs < 0) {
    fprintf(errfp, MSG_KINS_NEG_MAXRS);
    return(LIN_ILL_INPUT);
  }

  kinspgmr_mem->g_maxlrst = maxrs;

  return(SUCCESS);
}


/**********************************************************************/

int KINSpgmrSetPrecSetupFn(void *kinmem, KINSpgmrPrecSetupFn pset)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  kinspgmr_mem->g_pset = pset;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrSetPrecSolveFn(void *kinmem, KINSpgmrPrecSolveFn psolve)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  kinspgmr_mem->g_psolve = psolve;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrSetPrecData(void *kinmem, void *P_data)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  kinspgmr_mem->g_P_data = P_data;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrSetJacTimesVecFn(void *kinmem, KINSpgmrJacTimesVecFn jtimes)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  kinspgmr_mem->g_jtimes = jtimes;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrSetJacData(void *kinmem, void *J_data)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  kinspgmr_mem->g_J_data = J_data;

  return(SUCCESS);
}

/**********************************************************************/

int KINSpgmrGetIntWorkSpace(void *kinmem, long int *leniwSG)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;
  int maxl;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  maxl = kinspgmr_mem->g_maxl;
  *leniwSG = liw1*(maxl+3);

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetRealWorkSpace(void *kinmem, long int *lenrwSG)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;
  int maxl;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  maxl = kinspgmr_mem->g_maxl;
  *lenrwSG = lrw1*(maxl+3) + maxl*(maxl+4) + 1;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumPrecEvals(void *kinmem, long int *npevals)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *npevals = npe;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumPrecSolves(void *kinmem, long int *npsolves)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *npsolves = nps;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumLinIters(void *kinmem, long int *nliters)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *nliters = nli;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumConvFails(void *kinmem, long int *nlcfails)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *nlcfails = ncfl;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumJtimesEvals(void *kinmem, long int *njvevals)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *njvevals = njtimes;

  return(OKAY);
}

/**********************************************************************/

int KINSpgmrGetNumFuncEvals(void *kinmem, long int *nfevalsSG)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_SETGET_KINMEM_NULL);
    return(LIN_NO_MEM);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  kinspgmr_mem = (KINSpgmrMem) lmem;

  *nfevalsSG = nfeSG;

  return(OKAY);
}


/* Additional readability Replacements */

#define maxl       (kinspgmr_mem->g_maxl)
#define maxlrst    (kinspgmr_mem->g_maxlrst)

#define pset       (kinspgmr_mem->g_pset)
#define psolve     (kinspgmr_mem->g_psolve)
#define P_data     (kinspgmr_mem->g_P_data)

#define jtimes     (kinspgmr_mem->g_jtimes)
#define J_data     (kinspgmr_mem->g_J_data)

/*************** KINSpgmrInit *****************************************

 This routine initializes variables associated with the GMRES linear
 solver. Memory allocation  was done previously in KINSPgmr.        

**********************************************************************/

static int KINSpgmrInit(KINMem kin_mem)
{
  KINSpgmrMem kinspgmr_mem;

  kinspgmr_mem = (KINSpgmrMem) lmem;


  /* Initialize counters */

  npe = nli = nps = ncfl = 0;
  njtimes = nfeSG = 0;

  /* Set preconditioner type */
  if (psolve != NULL) {
    pretype = RIGHT;
  } else {
    pretype = NONE;
  }
  
  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = KINSpgmrDQJtimes;
    J_data = kin_mem;
  }

  return(LINIT_OK);
}

/*************** KINSpgmrSetup ****************************************

 This routine does the setup operations for the Spgmr linear solver, that is,
 it is an interface to the user-supplied routine  precondset 
**********************************************************************/

static int KINSpgmrSetup(KINMem kin_mem)
{
  
  int  ret;
  KINSpgmrMem kinspgmr_mem;

  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Call pset routine */
  ret = pset(uu, uscale, fval, fscale, P_data, vtemp1, vtemp2); 

  if(ret != 0) return(1);
  npe++;
  nnilpre = nni; 
  precondcurrent = TRUE; 

  /* Return the same value ret that precondset returned */
  return(0);
}

/*************** KINSpgmrSolve ****************************************

 This routine handles the call to the generic SPGMR solver SpgmrSolve
 for the solution of the linear system Ax = b.

 Appropriate variables are passed to SpgmrSolve and the counters 
 nli, nps, and ncfl are incremented, and the return value is set 
 according to the success of SpgmrSolve.  The success flag is
 returned if SpgmrSolve converged, or if this is the first Newton
 iteration and the residual norm was reduced below its initial value.
 Of the other error conditions, only preconditioner solver failure is
 specifically returned.  Otherwise a generic flag is returned to
 denote failure of this routine.

**********************************************************************/

static int KINSpgmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb, 
                         realtype *res_norm)
{
  KINSpgmrMem kinspgmr_mem;
  int ret;
  int nli_inc, nps_inc;
  
  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Set input initial guess xx = 0 to SpgmrSolve.   bb is set, by the routine
     calling KINSpgmrSolve, to the RHS vector for the system to be solved*/ 
 
  N_VConst(ZERO, xx);

  new_uu = TRUE;  /* set flag required for user Jacobian rtn */

  /* Call SpgmrSolve  */
  ret = SpgmrSolve(spgmr_mem, kin_mem, xx, bb, pretype, gstype, eps, 
                   maxlrst, kin_mem, fscale, fscale, KINSpgmrAtimes,
                   KINSpgmrPSolve, res_norm, &nli_inc, &nps_inc);
  /* Increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */
  nli = nli + (long int) nli_inc;
  nps = nps + (long int) nps_inc;

  if (printfl == 3) fprintf(errfp,"KINSpgmrSolve: nli_inc=%d\n",nli_inc);

  if (ret != 0) ncfl++;

  /* Compute the terms sJpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm.  Both of these terms are subsequently
     corrected if the step is reduced by constraints or the linesearch.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale.                        */

  KINSpgmrAtimes(kin_mem, xx, bb);
  sJpnorm = N_VWL2Norm(bb,fscale);
  N_VProd(bb, fscale, bb);
  N_VProd(bb, fscale, bb);
  sfdotJp = N_VDotProd(fval, bb);

  if (printfl > TWO) fprintf(errfp,
     "linear (Krylov step) residual norm=%12.3g  eps=%12.3g\n",*res_norm, eps);

  /* Set return value to appropriate values */

 if ( ret == SPGMR_SUCCESS || ret == SPGMR_RES_REDUCED ) return(0);

 if ( ret == SPGMR_PSOLVE_FAIL_REC ) return(1);

 if (ret == SPGMR_PSOLVE_FAIL_UNREC)
   return(KINSOL_PRECONDSOLVE_FAILURE);
 else
   return(KINSOL_KRYLOV_FAILURE);  
}

/*************** KINSpgmrFree ****************************************

 This routine frees memory specific to the Spgmr linear solver.

**********************************************************************/

static int KINSpgmrFree(KINMem kin_mem)
{
  KINSpgmrMem kinspgmr_mem;

  kinspgmr_mem = (KINSpgmrMem) lmem;

  SpgmrFree(spgmr_mem);
  free(lmem);
  return(0);
}

/*************** KINSpgmrAtimes *************************************

 This routine coordinates the generation of the matrix-vector product 
                       z = J * v 
 by calling either KINSpgmrAtimesDQ, which uses a difference quotient
 approximation for J*v, or by calling the user supplied routine
 userAtimes if it is non-null.

**********************************************************************/

static int KINSpgmrAtimes(void *kinsol_mem, N_Vector v, N_Vector z)
{
  KINMem   kin_mem;
  KINSpgmrMem kinspgmr_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinspgmr_mem = (KINSpgmrMem) lmem;

  ret = jtimes(v, z, uu, &new_uu, J_data);
  njtimes++;

  return(ret);
}
/*************** KINSpgmrPSolve ***************************************

 This routine interfaces between the generic SpgmrSolve routine and
 the user's psolve routine.  It passes to psolve all required state 
 information from kinsol_mem.  Its return value is the same as that
 returned by psolve. Note that the generic SPGMR solver guarantees
 that KINSpgmrPSolve will not be called in the case in which
 preconditioning is not done. This is the only case in which the
 user's psolve routine is allowed to be NULL.

**********************************************************************/

static int KINSpgmrPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lrdummy)
{
  KINMem   kin_mem;
  KINSpgmrMem kinspgmr_mem;
  int ret;

  kin_mem = (KINMem) kinsol_mem;
  kinspgmr_mem = (KINSpgmrMem)lmem;

  /* Copy the rhs into z before the psolve call;   
              z returns with the solution */
  N_VScale(ONE, r, z);

  ret = psolve(uu, uscale, fval, fscale, z, P_data, vtemp1);
  
  /* This call is counted in nps within the KINSpgmrSolve routine */

  return(ret);     
}

/*************** KINSpgmrDQJtimes *************************************

 This routine generates the matrix-vector product 
                       z = J * v 
 by using a difference quotient approximation.  The approximation is 
 J * v = [func(uu + sigma*v) - func(uu)]/sigma.  Here sigma is based on
 the dot products (uscale*uu, uscale*v) and (uscale*v, uscale*v),
 the L1Norm(uscale*v), and on sqrt_relfunc (the square root of the
 relative error in the function).  Note that v in the argument list
 has already been both preconditioned and unscaled.

**********************************************************************/

static int KINSpgmrDQJtimes(N_Vector v, N_Vector Jv,
                            N_Vector u, booleantype *new_u, 
                            void *jac_data)
{
  KINMem   kin_mem;
  KINSpgmrMem kinspgmr_mem;
  realtype sigma, sigma_inv;
  realtype sutsv, sq1norm, sign, vtv;

  /* jac_data is kin_mem */

  kin_mem = (KINMem) jac_data;
  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Scale the vector v , put Du*v into vtemp1. */
  N_VProd(v, uscale, vtemp1);
  /*  Scale u and put into Jv, used as a temporary. */
  N_VProd(u, uscale, Jv);

  /* Compute dot product (Du*u).(Du*v). */
  sutsv = N_VDotProd(Jv, vtemp1);

  /* Compute dot product (Du*v).(Du*v). */
  vtv = N_VDotProd(vtemp1, vtemp1);

  sq1norm = N_VL1Norm(vtemp1);

  sign = (sutsv >= ZERO) ? ONE : -ONE ;
 
  /*  This expression for sigma is from p. 469, Brown and Saad paper. */
  sigma = sign*sqrt_relfunc*MAX(ABS(sutsv),sq1norm)/vtv; 

  sigma_inv = ONE/sigma;

  /*  Compute the u-prime at which to evaluate the function func. */
  N_VLinearSum(ONE, u, sigma, v, vtemp1);
 
  /*  Call the system function to calc func(u+sigma*v).  */
  func(vtemp1, vtemp2, f_data);    
  nfeSG++;

  /*  Finish the computation of the difference quotient. */
  N_VLinearSum(sigma_inv, vtemp2, -sigma_inv, fval, Jv);

  return(0);
}
