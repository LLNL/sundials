/*******************************************************************
 *                                                                 *
 * File          : kinspgmr.c                                      *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 26 July 2002                                    *
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
#include "iterativ.h"
#include "spgmr.h"


/* Error Messages */

#define KINSPGMR_MAIN  "KINSpgmr-- "

#define MSG_MEM_FAIL    KINSPGMR_MAIN "A memory request failed.\n\n"

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/******************************************************************
 *                                                                *           
 * Types : KINSpgmrMemRec, KINSpgmrMem                            *
 *----------------------------------------------------------------*
 * The type KINSpgmrMem is pointer to a KINSpgmrMemRec. This      *
 * structure contains KINSpgmr solver-specific data.              *
 *                                                                *
 ******************************************************************/

typedef struct {

  int  g_maxl;        /* maxl = maximum dimension of the Krylov space   */     
  int  g_pretype;     /* preconditioning type--for Spgmr call           */
  int  g_gstype;      /* gram schmidt type --  for Spgmr call           */
  booleantype g_new_uu;  /* flag that a new uu has been created--
                            indicating that a call to generate a new
                            user-supplied Jacobian is required */
  int g_maxlrst;      /* max number of linear solver restarts allowed;
                         default is zero  */
  long int g_nli;     /* nli = total number of linear iterations        */
  long int g_npe;     /* npe = total number of precondset calls         */
  long int g_nps;     /* nps = total number of psolve calls             */
  long int g_ncfl;    /* ncfl = total number of convergence failures    */

  KINSpgmrPrecondFn g_precond; 
                      /* precondset = user-supplied routine to
                         compute a preconditioner                       */

  KINSpgmrPrecondSolveFn g_psolve; 
                      /* psolve = user-supplied routine to 
                         solve preconditioner linear system             */ 

  KINSpgmruserAtimesFn g_userAtimes;
                      /* userAtimes = user-supplied routine to optionally
                         compute the product J v as required by Spgmr   */ 

  void *g_P_data;     /* P_data is a memory block passed to psolve 
                         and precondset through KINSpgmr                */
  SpgmrMem g_spgmr_mem;  
                      /* spgmr_mem is memory used by the        
                         generic Spgmr solver                           */

} KINSpgmrMemRec, *KINSpgmrMem;


/* KINSpgmr linit, lsetup, lsolve, and lfree routines */

static int KINSpgmrInit(KINMem kin_mem, booleantype *setupNonNull);

static int KINSpgmrSetup(KINMem kin_mem);

static int KINSpgmrSolve(KINMem kin_mem, N_Vector xx, N_Vector bb,
                         realtype *res_norm);

static int KINSpgmrFree(KINMem kin_mem);

/* KINSpgmr Atimes and PSolve routines called by generic SPGMR solver */

static int KINSpgmrAtimes(void *kinsol_mem, N_Vector v, N_Vector z);

static int KINSpgmrAtimesDQ(void *kinsol_mem, N_Vector v, N_Vector z);

static int KINSpgmrPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lr);


/* Readability Replacements */

#define Neq     (kin_mem->kin_Neq)      
#define uround  (kin_mem->kin_uround)
#define nfe     (kin_mem->kin_nfe)
#define nni     (kin_mem->kin_nni)
#define nnilpre (kin_mem->kin_nnilpre)
#define func    (kin_mem->kin_func)
#define f_data  (kin_mem->kin_f_data)
#define msgfp   (kin_mem->kin_msgfp)
#define printfl (kin_mem->kin_printfl)
#define iopt    (kin_mem->kin_iopt)
#define linit   (kin_mem->kin_linit)
#define lsetup  (kin_mem->kin_lsetup)
#define lsolve  (kin_mem->kin_lsolve)
#define lfree   (kin_mem->kin_lfree)
#define lmem    (kin_mem->kin_lmem)
#define machenv (kin_mem->kin_machenv)
#define uu      (kin_mem->kin_uu)
#define fval    (kin_mem->kin_fval)
#define uscale  (kin_mem->kin_uscale)
#define fscale  (kin_mem->kin_fscale)
#define sqrt_relfunc (kin_mem->kin_sqrt_relfunc)
#define precondflag (kin_mem->kin_precondflag)

#define spgmr_mem (kinspgmr_mem->g_spgmr_mem)
#define nli     (kinspgmr_mem->g_nli)
#define npe     (kinspgmr_mem->g_npe)
#define nps     (kinspgmr_mem->g_nps)
#define ncfl    (kinspgmr_mem->g_ncfl)


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

int KINSpgmr(void *kinsol_mem, int maxl, int maxlrst, int msbpre,
             KINSpgmrPrecondFn precondset, KINSpgmrPrecondSolveFn psolve,
             KINSpgmruserAtimesFn userAtimes, void *P_data)
{
  KINMem kin_mem;
  KINSpgmrMem kinspgmr_mem;

  kin_mem = (KINMem)kinsol_mem;

  if (kin_mem == NULL){
     fprintf(msgfp, MSG_MEM_FAIL);
     return(KIN_MEM_NULL);  
  }

  /* Set four main function fields in kin_mem */
  linit  = KINSpgmrInit; 
  lsetup = KINSpgmrSetup;
  lsolve = KINSpgmrSolve;
  lfree  = KINSpgmrFree;


  /* Get memory for KINSpgmrMemRec */
  lmem = kinspgmr_mem = (KINSpgmrMem) malloc(sizeof(KINSpgmrMemRec));
  if (kinspgmr_mem == NULL){
     fprintf(msgfp, MSG_MEM_FAIL);
     return(KINSPGMR_MEM_FAIL);  
  }

  /*  Set Spgmr parameters appropriately for this package. The only choices
      with repect to preconditioning are NONE or RIGHT. The other options are 
      not available as they were with CVODE/PVODE, where pretype and gstype 
      were inputs to CVSpgmr. Here, the choice of NONE or RIGHT is implied 
      by the  'NULL'ness of the pointer to the psolve routine.              */

  if(psolve == NULL){
    precondflag = FALSE;
    kinspgmr_mem->g_pretype = NONE;
    }
  else{
    precondflag = TRUE;
    kinspgmr_mem->g_pretype = RIGHT;
    }

  kinspgmr_mem->g_gstype  = MODIFIED_GS;

  /* Set Spgmr parameters that have been passed in call sequence */
  kinspgmr_mem->g_maxl = (maxl<= 0) ? MIN(KINSPGMR_MAXL,Neq) : MIN(maxl,Neq);
  kinspgmr_mem->g_maxlrst = (maxlrst<=0) ? 0 : 
                             MIN(maxlrst,2*Neq/(kinspgmr_mem->g_maxl));
  kinspgmr_mem->g_P_data  = P_data;
  kinspgmr_mem->g_precond = precondset;
  kinspgmr_mem->g_psolve  = psolve;
  kinspgmr_mem->g_userAtimes = userAtimes;


  /*  This sets a variable in the main memory block for use in controlling
      the calling of the preconditioner. It is set to either the default 
      KINSPGMR_SBPRE or the user supplied value msbpre.                  */
  kin_mem->kin_msbpre  = (msbpre<=0) ? KINSPGMR_MSBPRE : msbpre ;

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(Neq, kinspgmr_mem->g_maxl, machenv);

  if (spgmr_mem == NULL) {
    fprintf(msgfp, MSG_MEM_FAIL);
    lmem = NULL; /* Set lmem to NULL and free that memory, as a flag to a
                    later inadvertent KINSol call, that SpgmrMalloc failed. */
    free(lmem);
    return(SPGMR_MEM_FAIL);
  }


  return(0);
}


/* Additional readability Replacements */
#define pretype    (kinspgmr_mem->g_pretype)
#define gstype     (kinspgmr_mem->g_gstype)
#define psolve     (kinspgmr_mem->g_psolve)
#define precondset (kinspgmr_mem->g_precond)
#define P_data     (kinspgmr_mem->g_P_data)
#define userAtimes (kinspgmr_mem->g_userAtimes)
#define maxl       (kinspgmr_mem->g_maxl)

/*************** KINSpgmrInit *****************************************

 This routine initializes variables associated with the GMRES linear
 solver. Memory allocation  was done previously in KINSPgmr.        

**********************************************************************/

static int KINSpgmrInit(KINMem kin_mem, booleantype *setupNonNull)
{
  KINSpgmrMem kinspgmr_mem;

  kinspgmr_mem = (KINSpgmrMem) lmem;


  /* Initialize counters, and set workspace lengths */

  npe = nli = nps = ncfl = 0;
    

  if (iopt != NULL) {
    iopt[SPGMR_NPE] = npe;
    iopt[SPGMR_NLI] = nli;
    iopt[SPGMR_NPS] = nps;
    iopt[SPGMR_NCFL] = ncfl;
  }


  /* Set setupNonNull to TRUE iff there is preconditioning:       */
  /* (pretype != NONE) and there is a preconditioning setup phase */
  /* (precond != NULL)                                            */

  *setupNonNull = ( (kinspgmr_mem->g_pretype  != NONE) && 
                    (kinspgmr_mem->g_precond  != NULL ) );

  return(LINIT_OK);
}

#define vtemp1  (kin_mem->kin_vtemp1)
#define vtemp2  (kin_mem->kin_vtemp2)
#define precondcurrent (kin_mem->kin_precondcurrent)

/*************** KINSpgmrSetup ****************************************

 This routine does the setup operations for the Spgmr linear solver, that is,
 it is an interface to the user-supplied routine  precondset 
**********************************************************************/

static int KINSpgmrSetup(KINMem kin_mem)
{
  
  int  ret;
  KINSpgmrMem kinspgmr_mem;

  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Call precondset routine */
  ret = precondset(Neq, uu, uscale, fval, fscale, vtemp1, vtemp2, func, 
        uround, &nfe, P_data);

  if(ret != 0) return(1);
  npe++;
  nnilpre = nni; 
  precondcurrent = TRUE; 

  /* Set npe, and return the same value ret that precondset returned */
  if (iopt != NULL) iopt[SPGMR_NPE] = npe;
  return(0);
}


/*   More readability constants defined  */

#define eps            (kin_mem->kin_eps)
#define ioptExists     (kin_mem->kin_ioptExists)
#define maxlinrestarts (kinspgmr_mem->g_maxlrst)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)


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
  int ret, nli_inc, nps_inc;
  
  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Set input initial guess xx = 0 to SpgmrSolve.   bb is set, by the routine
     calling KINSpgmrSolve, to the RHS vector for the system to be solved*/ 
 
  N_VConst(ZERO, xx);

  kinspgmr_mem->g_new_uu = TRUE;  /* set flag required for user Jacobian rtn */

  /* Call SpgmrSolve  */
  ret = SpgmrSolve(spgmr_mem, kin_mem, xx, bb, pretype, gstype, eps, 
                   maxlinrestarts, kin_mem, fscale, fscale, KINSpgmrAtimes,
                   KINSpgmrPSolve, res_norm, &nli_inc, &nps_inc);
  /* Increment counters nli, nps, and ncfl 
     (nni is updated in the KINSol main iteration loop) */
  nli += nli_inc;
  nps += nps_inc;

  if (printfl == 3) fprintf(msgfp,"KINSpgmrSolve: nli_inc=%d\n",nli_inc);

  if (ioptExists) {
    iopt[SPGMR_NLI] = nli;
    iopt[SPGMR_NPS] = nps;
  }  
  if (ret != 0) { 
    ncfl++;
    if (ioptExists) iopt[SPGMR_NCFL] = ncfl;
  }

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

  if (printfl > TWO) fprintf(msgfp,
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
  integertype ret;
  KINSpgmruserAtimesFn patimes;

  kin_mem = (KINMem) kinsol_mem;
  kinspgmr_mem = (KINSpgmrMem) lmem;

  patimes = kinspgmr_mem->g_userAtimes;

  if(patimes == NULL)
    ret = KINSpgmrAtimesDQ(kinsol_mem, v,z);
  else 
    ret = (*userAtimes)((kin_mem->kin_f_data), v, z, 
                        &(kinspgmr_mem->g_new_uu),uu);
  
  return(ret);
}
/*************** KINSpgmrAtimesDQ *************************************

 This routine generates the matrix-vector product 
                       z = J * v 
 by using a difference quotient approximation.  The approximation is 
 J * v = [func(uu + sigma*v) - func(uu)]/sigma.  Here sigma is based on
 the dot products (uscale*uu, uscale*v) and (uscale*v, uscale*v),
 the L1Norm(uscale*v), and on sqrt_relfunc (the square root of the
 relative error in the function).  Note that v in the argument list
 has already been both preconditioned and unscaled.

**********************************************************************/

static int KINSpgmrAtimesDQ(void *kinsol_mem, N_Vector v, N_Vector z)
{
  realtype sigma, sigma_inv;
  realtype sutsv, sq1norm, sign, vtv;
  KINMem   kin_mem;
  KINSpgmrMem kinspgmr_mem;

  kin_mem = (KINMem) kinsol_mem;
  kinspgmr_mem = (KINSpgmrMem) lmem;

  /* Scale the vector v , put Du*v into vtemp1. */
  N_VProd(v, uscale, vtemp1);
  /*  Scale uu and put into z, used as a temporary. */
  N_VProd(uu, uscale, z);

  /* Compute dot product (Du*u).(Du*v). */
  sutsv = N_VDotProd(z, vtemp1);

  /* Compute dot product (Du*v).(Du*v). */
  vtv = N_VDotProd(vtemp1, vtemp1);

  sq1norm = N_VL1Norm(vtemp1);

  sign = (sutsv >= ZERO) ? ONE : -ONE ;
 
  /*  This expression for sigma is from p. 469, Brown and Saad paper. */
  sigma = sign*sqrt_relfunc*MAX(ABS(sutsv),sq1norm)/vtv; 

  sigma_inv = ONE/sigma;

  /*  Compute the u-prime at which to evaluate the function func. */
  N_VLinearSum(ONE, uu, sigma, v, vtemp1);
 
  /*  Call the system function to calc func(uu+sigma*v).  */
  func(Neq, vtemp1, vtemp2, f_data);    nfe++;

  /*  Finish the computation of the difference quotient. */
  N_VLinearSum(sigma_inv, vtemp2, -sigma_inv, fval, z);

  return(0);
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

  ret = psolve(Neq, uu, uscale, fval, fscale, z, vtemp1, func, uround,
               &nfe, P_data);
  
  /* This call is counted in nps within the KINSpgmrSolve routine */

  return(ret);     
}
