/*******************************************************************
 * File          : kinsol.c                                        *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/kinsol/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for the main KINSol solver.     *
 * It is independent of the KINSol linear solver in use.           *
 *                                                                 *
 *******************************************************************/

/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "kinsol.h"
#include "kinspgmr.h"

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/

/***************************************************************/
/*********************** BEGIN Macros **************************/
/***************************************************************/
/* Macro: loop */

#define loop for(;;)
/***************************************************************/
/************************ END Macros ***************************/
/***************************************************************/

/************************************************************/
/************** BEGIN KINSol Private Constants ***************/
/************************************************************/

#define HALF   RCONST(0.5)      /* real 0.5   */
#define ZERO   RCONST(0.0)      /* real 0.0   */
#define ONE    RCONST(1.0)      /* real 1.0   */
#define ONEPT5 RCONST(1.5)      /* real 1.5   */
#define TWO    RCONST(2.0)      /* real 2.0   */
#define FIVE   RCONST(5.0)      /* real 5.0   */
#define TEN    RCONST(10.0)     /* real 10.0  */
#define POINT1 RCONST(0.1)      /* real 0.1   */
#define POINTOH1 RCONST(0.01)   /* real 0.01  */
#define POINT99  RCONST(0.99)   /* real 0.99  */
#define THOUSAND RCONST(1000.0) /* real 1000.0 */
#define ONETHIRD RCONST(.3333333333333333)    /* 1/3  */
#define TWOTHIRDS RCONST(.6666666666666667)   /* 2/3  */
#define POINT6   RCONST(0.6)    /* real 0.6  */
#define POINT9   RCONST(0.9)    /* real 0.9  */
#define POINTOHOHONE RCONST(0.001)    /* real 0.001 */
#define POINTOHOHOHONE RCONST(0.0001) /* real 0.0001 */

/***************************************************************/
/************** BEGIN Default Constants ************************/
/***************************************************************/

#define RELU_DEFAULT RCONST(10000.0) 
#define PRINTFL_DEFAULT 0
#define MXITER_DEFAULT 200
#define MXNBCF 10
#define MSBPRE 10

/***************************************************************/
/*************** END Default Constants *************************/
/***************************************************************/

/* KINStop return value requesting more iterations */
#define CONTINUE_ITERATIONS -999

/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* KINCreate Error Messages */

#define MSG_KINMEM_FAIL   "KINCreate-- Allocation of kin_mem failed.\n\n "

/* KINSet* Error Messages */

#define MSG_KINS_NO_MEM   "kin_mem=NULL in a KINSet routine illegal.\n\n"

#define MSG_BAD_PRINTFL   "KINSetPrintLevel-- illegal value for printfl.\n\n"

#define MSG_BAD_MXITER    "KINSetNumMaxIters-- illegal value for mxiter.\n\n"

#define MSG_BAD_MSBPRE    "KINSetMaxPrecCalls-- illegal msbpre<0. \n\n"

#define MSG_BAD_ETACHOICE "KINSetEtaForm-- illegal value for etachoice.\n\n"

#define MSG_BAD_ETACONST  "KINSetEtaConstValue-- eta out of range.\n\n"

#define MSG_BAD_GAMMA     "KINSetEtaParams-- gamma out of range.\n\n"
#define MSG_BAD_ALPHA     "KINSetEtaParams-- alpha out of range.\n\n"

#define MSG_BAD_MXNEWTSTEP "KINSetMaxNewtonStep-- mxnewtstep nonpositive.\n\n"

#define MSG_BAD_RELFUNC   "KINSetRelErrFunc-- relfunc=%g < 0 illegal.\n\n"

#define MSG_BAD_RELU      "KINSetMaxSolUpdate-- relu=%g < 0 illegal.\n\n"

#define MSG_BAD_FNORMTOL  "KINSetFuncNormTol-- fnormtol=%g < 0 illegal.\n\n"
#define MSG_BAD_SCSTEPTOL "KINSetScaledStepTol-- scsteptol=%g < 0 illegal.\n\n"

/* KINMalloc Error Messages */

#define MSG_KINM_NO_MEM   "KINMalloc-- kin_mem=NULL illegal.\n\n"

#define MSG_MEM_FAIL      "KINMalloc-- A memory request failed.\n\n"

#define MSG_FUNC_NULL     "KINMalloc-- func=NULL illegal.\n\n"

#define MSG_KINS_FUNC_NULL "KINResetSysFunc-- func=NULL illegal.\n\n"

/* KINSol error messages */

#define MSG_KINSOL_NO_MEM    "KINSol-- kinsol_mem=NULL illegal.\n\n"

#define MSG_KINSOL_NO_MALLOC "KINSol-- Attempt to call before KINMalloc illegal.\n\n"

/* KINSolInit messages */

#define KINSI "KINSolInit--"

#define MSG_LSOLV_NO_MEM    KINSI "The linear solver memory pointer is NULL.\n\n"
 
#define MSG_UU_NULL         KINSI "uu=NULL illegal.\n\n"

#define MSG_BAD_GLSTRAT1    KINSI "globalstrategy=%d illegal.\n"
#define MSG_BAD_GLSTRAT2    "The legal values are INEXACT_NEWTON=%d"
#define MSG_BAD_GLSTRAT3    " and LINESEARCH=%d\n\n"
#define MSG_BAD_GLSTRAT     MSG_BAD_GLSTRAT1 MSG_BAD_GLSTRAT2 MSG_BAD_GLSTRAT3

#define MSG_BAD_USCALE      KINSI "uscale=NULL illegal.\n\n"
#define MSG_USCALE_NONPOSITIVE  KINSI "uscale has nonpositive elements.\n\n"

#define MSG_BAD_FSCALE      KINSI "fscale=NULL illegal.\n\n"
#define MSG_FSCALE_NONPOSITIVE  KINSI "fscale has nonpositive elements.\n\n"
 
#define MSG_INITIAL_CNSTRNT KINSI "Initial guess does NOT meet constraints.\n\n"

#define MSG_LINIT_FAIL      KINSI "The linear solver's init routine failed.\n\n"

/* KINGet* Error Messages */

#define MSG_KING_NO_MEM     "kin_mem=NULL in a KINGGet* routine illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/

/**************************************************************/
/********* BEGIN Private Helper Functions Prototypes **********/
/**************************************************************/

static booleantype KINAllocVectors(KINMem kin_mem, NV_Spec NVSpec);

static int KINSolInit(KINMem kin_mem);

static int KINConstraint(KINMem kin_mem );

static void KINForcingTerm(KINMem kin_mem, realtype fnormp);

static void KINFreeVectors(KINMem kin_mem);

static int  KINInexactNewton(KINMem kin_mem, realtype *fnormp, 
                             realtype *f1normp, booleantype *maxStepTaken);

static int  KINLineSearch(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxSteptaken);

static booleantype KINInitialStop(KINMem kin_mem);

static int  KINLinSolDrv(KINMem kinmem , N_Vector bb , N_Vector xx );

static realtype KINScFNorm(N_Vector vv , N_Vector scale , N_Vector wrkv);

static realtype KINScSteplength(KINMem kin_mem,
                                N_Vector ucur, N_Vector ss, N_Vector usc);

static int KINStop(KINMem kinmem, booleantype maxStepTaken, int globalstratret);

/**************************************************************/
/********** END Private Helper Functions Prototypes ***********/
/**************************************************************/

/***************************************************************/
/************* BEGIN KINSol Implementation **********************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/

/*------------------       KINCreate     --------------------------*/
/* 
   KINCreate creates an internal memory block for a problem to 
   be solved by KINSOL.
   If successful, KINCreate returns a pointer to the problem memory. 
   This pointer should be passed to KINMalloc.  
   If an initialization error occurs, KINCreate prints an error 
   message to standard err and returns NULL. 
*/
/*-----------------------------------------------------------------*/

void *KINCreate(void)
{
  KINMem kin_mem;
  realtype uround;

  kin_mem = (KINMem) malloc(sizeof(struct KINMemRec));
  if (kin_mem == NULL) {
    fprintf(stdout, MSG_KINMEM_FAIL);
    return(NULL);
  }

  /* Set uround */
  kin_mem->kin_uround = uround = UNIT_ROUNDOFF;
  
  /* Set default values for solver optional inputs */
  kin_mem->kin_f_data = NULL;
  kin_mem->kin_errfp = stdout;
  kin_mem->kin_printfl = PRINTFL_DEFAULT;
  kin_mem->kin_mxiter = MXITER_DEFAULT;
  kin_mem->kin_noPrecInit = FALSE;
  kin_mem->kin_msbpre = MSBPRE;
  kin_mem->kin_pthrsh = TWO;
  kin_mem->kin_noMinEps = FALSE;
  kin_mem->kin_mxnewtstep = ZERO;
  kin_mem->kin_sqrt_relfunc = RSqrt(uround);
  kin_mem->kin_relu = RELU_DEFAULT;
  kin_mem->kin_scsteptol = RPowerR(uround,TWOTHIRDS);
  kin_mem->kin_fnormtol = RPowerR(uround,ONETHIRD);

  kin_mem->kin_etaflag   = ETACHOICE1;
  kin_mem->kin_eta       = POINT1;     /* default for ETACONSTANT */
  kin_mem->kin_eta_alpha = TWO;        /* default for ETACHOICE2 */
  kin_mem->kin_eta_gamma = POINT9;     /* default for ETACHOICE2 */


  kin_mem->kin_MallocDone = FALSE;

  return((void *)kin_mem);
}

/*=================================================================*/
/*BEGIN        SOLVER OPTIONAL INPUT FUNCTIONS                     */
/*=================================================================*/

/*-----------------------------------------------------------------*/

int KINSetFdata(void *kinmem, void *f_data)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_f_data = f_data;

  return(SUCCESS);
}

#define f_data (kin_mem->kin_f_data)

/*-----------------------------------------------------------------*/

int KINSetErrFile(void *kinmem, FILE *errfp)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_errfp = errfp;

  return(SUCCESS);
}

#define errfp (kin_mem->kin_errfp)

/*-----------------------------------------------------------------*/

int KINSetPrintLevel(void *kinmem, int printfl)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (printfl<0 || printfl>3) {
    fprintf(errfp, MSG_BAD_PRINTFL);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_printfl = printfl;

  return(SUCCESS);
}

#define printfl (kin_mem->kin_printfl)

/*-----------------------------------------------------------------*/

int KINSetNumMaxIters(void *kinmem, int mxiter)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (mxiter <= 0) {
    fprintf(errfp, MSG_BAD_MXITER);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_mxiter = mxiter;

  return(SUCCESS);
}

#define mxiter (kin_mem->kin_mxiter)

/*-----------------------------------------------------------------*/

int KINSetNoPrecInit(void *kinmem, booleantype noPrecInit)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_noPrecInit = noPrecInit;

  return(SUCCESS);
}

#define noPrecInit (kin_mem->kin_noPrecInit)

/*-----------------------------------------------------------------*/

int KINSetMaxPrecCalls(void *kinmem, int msbpre)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (msbpre < 0) {
    fprintf(errfp, MSG_BAD_MSBPRE);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_msbpre = msbpre;

  return(SUCCESS);
}

#define msbpre (kin_mem->kin_msbpre)

/*-----------------------------------------------------------------*/

int KINSetEtaForm(void *kinmem, int etachoice)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if ( (etachoice != ETACONSTANT) && 
       (etachoice != ETACHOICE1)  && 
       (etachoice != ETACHOICE2) ) {
    fprintf(errfp, MSG_BAD_ETACHOICE);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_etaflag = etachoice;

  return(SUCCESS);
}

#define etaflag (kin_mem->kin_etaflag)

/*-----------------------------------------------------------------*/

int KINSetEtaConstValue(void *kinmem, realtype eta)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if ( (eta <= ZERO) || (eta > ONE) ) {
    fprintf(errfp, MSG_BAD_ETACONST);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_eta = eta;

  return(SUCCESS);
}

#define eta (kin_mem->kin_eta)

/*-----------------------------------------------------------------*/

int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (ealpha <= ONE || ealpha > TWO)
    if (ealpha != ZERO) {
      fprintf(errfp,MSG_BAD_ALPHA);
      return(KINS_ILL_INPUT);
    }
  
  if (ealpha == ZERO) 
    kin_mem->kin_eta_alpha = TWO;
  else
    kin_mem->kin_eta_alpha = ealpha;

  if (egamma <= ZERO || egamma > ONE)
    if (egamma != ZERO) {
      fprintf(errfp,MSG_BAD_GAMMA);
      return(KINS_ILL_INPUT);
    }

  if (egamma == ZERO)
    kin_mem->kin_eta_gamma = POINT9;
  else
    kin_mem->kin_eta_gamma = egamma;

  return(SUCCESS);
}

#define ealpha (kin_mem->kin_eta_alpha)
#define egamma (kin_mem->kin_eta_gamma)

/*-----------------------------------------------------------------*/

int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_noMinEps = noMinEps;

  return(SUCCESS);
}

#define noMinEps (kin_mem->kin_noMinEps)

/*-----------------------------------------------------------------*/

int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (mxnewtstep <= ZERO) {
    fprintf(errfp, MSG_BAD_MXNEWTSTEP);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_mxnewtstep = mxnewtstep;

  return(SUCCESS);
}

#define mxnewtstep (kin_mem->kin_mxnewtstep)

/*-----------------------------------------------------------------*/

int KINSetRelErrFunc(void *kinmem, realtype relfunc)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (relfunc <= ZERO) {
    fprintf(errfp, MSG_BAD_RELFUNC, relfunc);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_sqrt_relfunc = RSqrt(relfunc);

  return(SUCCESS);
}

#define relfunc (kin_mem->kin_sqrt_relfunc)

/*-----------------------------------------------------------------*/

int KINSetMaxSolUpdate(void *kinmem, realtype relu)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (relu <= ZERO) {
    fprintf(errfp, MSG_BAD_RELU, relu);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_relu = relu;

  return(SUCCESS);
}

#define relu (kin_mem->kin_relu)

/*-----------------------------------------------------------------*/

int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (fnormtol <= ZERO) {
    fprintf(errfp, MSG_BAD_FNORMTOL, fnormtol);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_fnormtol = fnormtol;

  return(SUCCESS);
}

#define fnormtol (kin_mem->kin_fnormtol)

/*-----------------------------------------------------------------*/

int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (scsteptol <= ZERO) {
    fprintf(errfp, MSG_BAD_SCSTEPTOL, scsteptol);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_scsteptol = scsteptol;

  return(SUCCESS);
}

#define scsteptol (kin_mem->kin_scsteptol)

/*-----------------------------------------------------------------*/

int KINSetConstraints(void *kinmem, N_Vector constraints)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  kin_mem->kin_constraints = constraints;

  return(SUCCESS);
}

#define constraints (kin_mem->kin_constraints)

/*=================================================================*/
/*END        SOLVER OPTIONAL INPUT FUNCTIONS                       */
/*=================================================================*/

/******************** KINMalloc************************************

 KINMalloc allocates memory for a problem or execution of KINSol. 
 If memory is successfully allocated, SUCCESS is returned.
 Otherwise, an error message is printed and an error flag returned.
 
*******************************************************************/

int KINMalloc(void *kinmem, SysFn func, NV_Spec nvspec)
{
  KINMem kin_mem;
  booleantype allocOK;
  long int liw1, lrw1;
  
  /* Check kinmem */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_KINM_NO_MEM);
    return(KINM_NO_MEM);
  }
  
  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    fprintf(errfp, MSG_FUNC_NULL);
    return(KINM_ILL_INPUT);
  }

  /* Set space requirements for one N_Vector */
  N_VSpace(nvspec, &lrw1, &liw1);
  kin_mem->kin_lrw1 = lrw1;
  kin_mem->kin_liw1 = liw1;  
  
  /* Allocate the vectors. */
  allocOK = KINAllocVectors(kin_mem, nvspec);
  if (!allocOK) {
    fprintf(errfp, MSG_MEM_FAIL);
    free(kin_mem);
    return(KINM_MEM_FAIL);
  }

  /* Copy the input parameters into KINSOL state */
  kin_mem->kin_func = func;
  kin_mem->kin_nvspec = nvspec;

  /* Set the linear solver addresses to NULL. */
  kin_mem->kin_linit = NULL;
  kin_mem->kin_lsetup = NULL;
  kin_mem->kin_lsolve = NULL;
  kin_mem->kin_lfree = NULL;
  kin_mem->kin_lmem = NULL;
  
  /* Problem memory has been successfully allocated. */
  kin_mem->kin_MallocDone = TRUE;
  return(SUCCESS);
}

#define nvspec (kin_mem->kin_nvspec)

/*-----------------------------------------------------------------*/

int KINResetSysFunc(void *kinmem, SysFn func)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KINS_NO_MEM);
    return(KINS_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    fprintf(errfp, MSG_KINS_FUNC_NULL);
    return(KINS_ILL_INPUT);
  }

  kin_mem->kin_func = func;

  return(SUCCESS);
}

#define func (kin_mem->kin_func)

/**************************************************************/
/*********** BEGIN More Readability Constants *****************/
/**************************************************************/

#define uround   (kin_mem->kin_uround)
#define nni      (kin_mem->kin_nni)
#define nfe      (kin_mem->kin_nfe)
#define nbcf     (kin_mem->kin_nbcf)  
#define nbktrk   (kin_mem->kin_nbktrk)
#define ncscmx   (kin_mem->kin_ncscmx)
#define stepl    (kin_mem->kin_stepl)
#define pthrsh   (kin_mem->kin_pthrsh)
#define linit    (kin_mem->kin_linit)
#define lsetup   (kin_mem->kin_lsetup)
#define lsolve   (kin_mem->kin_lsolve) 
#define lfree    (kin_mem->kin_lfree)
#define constraintsSet (kin_mem->kin_constraintsSet) 
#define precondcurrent (kin_mem->kin_precondcurrent)          
#define nnilpre  (kin_mem->kin_nnilpre)
#define lmem     (kin_mem->kin_lmem)        
#define setupNonNull (kin_mem->kin_setupNonNull)
#define fval     (kin_mem->kin_fval)      
#define fnorm    (kin_mem->kin_fnorm)
#define f1norm   (kin_mem->kin_f1norm)
#define etaflag  (kin_mem->kin_etaflag)
#define callForcingTerm (kin_mem->kin_callForcingTerm)
#define uu        (kin_mem->kin_uu)
#define uscale    (kin_mem->kin_uscale)
#define fscale    (kin_mem->kin_fscale)
#define globalstrategy (kin_mem->kin_globalstrategy)     
#define sJpnorm   (kin_mem->kin_sJpnorm)
#define sfdotJp   (kin_mem->kin_sfdotJp)
#define unew      (kin_mem->kin_unew)
#define pp        (kin_mem->kin_pp)
#define vtemp1    (kin_mem->kin_vtemp1)
#define vtemp2    (kin_mem->kin_vtemp2)
#define eps       (kin_mem->kin_eps)
#define res_norm  (kin_mem->kin_res_norm)
#define precondflag (kin_mem->kin_precondflag)
#define liw1      (kin_mem->kin_liw1)
#define lrw1      (kin_mem->kin_lrw1)
#define liw       (kin_mem->kin_liw)
#define lrw       (kin_mem->kin_lrw)

/**************************************************************/
/************ END More Readability Constants ******************/
/**************************************************************/


/********************* KINSol *********************************************

   This routine is the main driver of the KINSol package. It manages the 
   computational process of computing an approximate solution to the
   system  func(uu) = 0. Routine KINLinSolDrv handles the process of 
   obtaining a solution to the system J x = b , where x is an increment
   from the last iterate uu. It then calls a global strategy routine 
   (either KINLinesearch or KINInexactNewton) to apply that increment xx 
   to the last iterate uu, creating a new iterate value uu. Finally 
   KINConstraint and KINStop are called to require that the user-specified 
   constraints are satisfied and if the iteration process has produced an 
   approximate iterate within the desired tolerance, resp.


***************************************************************************/

int KINSol(void *kinmem, N_Vector u, int strategy,  
           N_Vector u_scale, N_Vector f_scale)
{
  realtype fnormp, f1normp, epsmin;
  int  ret, globalstratret = 0;
  booleantype maxStepTaken;
  N_Vector bb, xx;
 
  KINMem kin_mem;

  /* Check for kinmem non NULL */
  if (kinmem == NULL) {
    fprintf(stdout, MSG_KINSOL_NO_MEM);
    return(KINSOL_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  if(kin_mem->kin_MallocDone == FALSE) {
    fprintf(errfp, MSG_KINSOL_NO_MALLOC);
    return(KINSOL_NO_MALLOC);
  }

  /* Load input arguments */
  globalstrategy = strategy;
  uu = u;
  uscale = u_scale;
  fscale = f_scale;

  /* Initialize solver */
  ret = KINSolInit(kin_mem);

  if (ret != SUCCESS) return(ret);
  
  ncscmx = 0;

  /* The following logic allows the choice of whether or not to force a
         call to the preconditioner precond upon a given call to KINSol. */

  if (noPrecInit == FALSE) pthrsh = TWO;

  if (noMinEps == FALSE) epsmin = POINTOH1*fnormtol;


  loop{

    nni++;

    /*  Calculate the epsilon for this iteration based on eta from the
           routine KINForcingTerm. */

    eps = (eta + uround) * fnorm;
    if(!noMinEps) eps = MAX(epsmin, eps);

    /*  Call KINLinSolDrv to calculate the approximate Newton step using
       the appropriate Krylov algorithm.  */

    /*  Use unew as work space for bb, the rhs vector, in the Newton steps 
        and pp as xx, the x in J x = b.  */

    bb = unew;  xx = pp;


    ret = KINLinSolDrv( kinmem, bb, xx );

    /* Process the return flag from KINLinSolDrv.  */

    if (ret != 0) break;

    /* Call the appropriate routine to calculate an acceptable step pp.  */

    if (globalstrategy == INEXACT_NEWTON) globalstratret =
      KINInexactNewton(kinmem, &fnormp, &f1normp, &maxStepTaken);

    else if (globalstrategy == LINESEARCH) globalstratret =
      KINLineSearch(kinmem, &fnormp, &f1normp, &maxStepTaken);
    
    /* If too many beta condition failures, stop iteration. */

    if (nbcf > MXNBCF) {
      ret = KINSOL_LINESEARCH_BCFAIL;
      break;
    }

    /* Evaluate eta by calling the forcing term routine. */

    if (callForcingTerm) KINForcingTerm(kinmem, fnormp);

    /*  Call KINStop to check if tolerances are met at this iteration. */

    fnorm = fnormp;

    ret = KINStop(kinmem, maxStepTaken, globalstratret); 

    /* Update uu after the iteration. */
    N_VScale(ONE, unew, uu);

    f1norm = f1normp;

    /*  Print out the current nni, fnorm, and nfe values if printfl > 0. */
    if (printfl>0) fprintf(errfp,"KINSol nni= %4d fnorm= %26.16g nfe=%6d\n",
                           nni, fnorm, nfe);

    if (ret != CONTINUE_ITERATIONS) break; 

    fflush(errfp);

 }   /* End of step loop; load optional outputs and return. */

  if(printfl>0){
    fprintf(errfp,"KINSol return value %d\n",ret);
    if(ret == SUCCESS)
        fprintf(errfp,"---SUCCESS\n");
    if(ret==KINSOL_STEP_LT_STPTOL)
        fprintf(errfp,"---KINSOL_STEP_LT_STPTOL\n");
    if(ret==KINSOL_LNSRCH_NONCONV)
        fprintf(errfp,"---KINSOL_LNSRCH_NONCONV");
    if(ret==KINSOL_LINESEARCH_BCFAIL)
        fprintf(errfp,"---KINSOL_LINESEARCH_BCFAIL");
    if(ret==KINSOL_MAXITER_REACHED)
        fprintf(errfp,"---KINSOL_MAXITER_REACHED\n");
    if(ret==KINSOL_MXNEWT_5X_EXCEEDED)
        fprintf(errfp,"---KINSOL_MXNEWT_5X_EXCEEDED\n");
    if(ret==KINSOL_KRYLOV_FAILURE)
        fprintf(errfp,"---KINSOL_KRYLOV_FAILURE\n");
    if(ret==KINSOL_PRECONDSET_FAILURE)
        fprintf(errfp,"---KINSOL_PRECONDSET_FAILURE\n");
    if(ret==KINSOL_PRECONDSOLVE_FAILURE)
        fprintf(errfp,"---KINSOL_PRECONDSOLVE_FAILURE\n");
  }

  return(ret);
}

/*-----------------------------------------------------------------*/

int KINGetIntWorkSpace(void *kinmem, long int *leniw)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *leniw = liw;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetRealWorkSpace(void *kinmem, long int *lenrw)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *lenrw = lrw;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetNumNonlinSolvIters(void *kinmem, int *nniters)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *nniters = nni;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetNumFuncEvals(void *kinmem, int *nfevals)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *nfevals = nfe;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetNumBetaCondFails(void *kinmem, int *nbcfails)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *nbcfails = nbcf;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetNumBacktrackOps(void *kinmem, int *nbacktr)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *nbacktr = nbktrk;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetFuncNorm(void *kinmem, realtype *funcnorm)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *funcnorm = fnorm;

  return(OKAY);
}

/*-----------------------------------------------------------------*/

int KINGetStepLength(void *kinmem, realtype *steplength)
{
  KINMem kin_mem;

  if (kinmem==NULL) {
    fprintf(stdout, MSG_KING_NO_MEM);
    return(KING_NO_MEM);
  }

  kin_mem = (KINMem) kinmem;

  *steplength = stepl;

  return(OKAY);
}

 
/********************* KINFree *************************************

  This routine frees the problem memory allocated by KINMalloc.
  Such memory includes all the vectors allocated by KINAllocVectors,
  and the memory lmem for the linear solver (deallocated by a call
  to lfree).

********************************************************************/

void KINFree(void *kinmem)
{
  KINMem kin_mem;

  if (kinmem == NULL) return;

  kin_mem = (KINMem) kinmem;

  KINFreeVectors(kin_mem);

  lfree(kin_mem);

  free(kin_mem);
}

/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Private Helper Functions Implementation ************/
/*******************************************************************/
 
/****************** KINAllocVectors ***********************************

 This routine allocates the KINSol vectors                           .
 The length of the vectors is the input
 parameter neq. If all memory allocations are successful,
 KINAllocVectors returns TRUE. Otherwise all allocated memory is freed
 and KINAllocVectors returns FALSE.

**********************************************************************/

static booleantype KINAllocVectors(KINMem kin_mem, NV_Spec NVSpec)
{
  /* Allocate unew, fval, pp, vtemp1 and vtemp2
     Any future modifier of this code is advised to be wary.
                  --- watch scope carefully  -- 
     unew, pp, vtemp1 and vtemp2 are used in  more than one context. */
  
  unew = N_VNew(NVSpec);
  if (unew == NULL) return(FALSE);

  fval = N_VNew(NVSpec);
  if (fval == NULL) {
    N_VFree(unew);
    return(FALSE);
  }

  pp = N_VNew(NVSpec);
  if (pp == NULL) {
    N_VFree(unew);
    N_VFree(fval);
    return(FALSE);
  }

  vtemp1 = N_VNew(NVSpec);
  if (vtemp1 == NULL) {
    N_VFree(unew);
    N_VFree(fval);
    N_VFree(pp);
    return(FALSE);
  }

  vtemp2 = N_VNew(NVSpec);
  if (vtemp2 == NULL) {
    N_VFree(unew);
    N_VFree(fval);
    N_VFree(pp);
    N_VFree(vtemp1);
    return(FALSE);
  }

  liw = 5*liw1;
  lrw = 5*lrw1;

  return(TRUE);
}

/******************** KINSolInit*********************************************

 KINSolInit initializes the problem for the specific input received in
 this call to KINSol (which calls KINSolInit). All problem specification 
 inputs are checked for errors. If any error occurs during initialization, 
 it is reported to the file whose file pointer is errfp.
 The possible return values for KINSolInit are:
    SUCCESS                 indicates a normal initialization. 
    KINSOL_LSOLV_NO_MEM     indicates that the linear solver has not been
                            specified
    KINSOL_INPUT_ERROR      indicates that an input error has been found
    KINSOL_INITIAL_GUESS_OK indicates that the guess uu satisfied the 
                            system func(uu) = 0 within the tolerances specified.

*****************************************************************************/

static int KINSolInit(KINMem kin_mem)
{
  int ret;
  
  /* Check for legal input parameters */

  if(lmem == NULL) {
    fprintf(errfp, MSG_LSOLV_NO_MEM);
    return(KINSOL_LSOLV_NO_MEM);
  }

  if (uu==NULL) {
    fprintf(errfp, MSG_UU_NULL);   
    return(KINSOL_INPUT_ERROR);
  }

  if ((globalstrategy != INEXACT_NEWTON) && (globalstrategy != LINESEARCH)) {
    fprintf(errfp, MSG_BAD_GLSTRAT, globalstrategy, INEXACT_NEWTON, LINESEARCH);
    return(KINSOL_INPUT_ERROR);
  }

  if (uscale == NULL)  {
    fprintf(errfp, MSG_BAD_USCALE);
    return(KINSOL_INPUT_ERROR);
  }

  if (N_VMin(uscale)<=ZERO){
    fprintf(errfp, MSG_USCALE_NONPOSITIVE); 
    return(KINSOL_INPUT_ERROR);
  }

  if (fscale == NULL)  {
    fprintf(errfp, MSG_BAD_FSCALE);
    return(KINSOL_INPUT_ERROR);
  }

  if (N_VMin(fscale)<=ZERO){
    fprintf(errfp, MSG_FSCALE_NONPOSITIVE);
    return(KINSOL_INPUT_ERROR);
  }

  /*  Set the constraints flag. */
  if (constraints == NULL) 
    constraintsSet = FALSE;
  else
    constraintsSet = TRUE;

  /* Check the initial guess uu against the constraints, if any. */
  if (constraintsSet) {
    if (!N_VConstrProdPos(constraints, uu)) {
      fprintf(errfp,MSG_INITIAL_CNSTRNT);
      KINFreeVectors(kin_mem);
      free(kin_mem);
      return(KINSOL_INPUT_ERROR);
    }
  }
  
  /* All error checking is complete at this point. */

  /* Copy the input parameters into KINSol state. */
  kin_mem->kin_uu = uu ;
  kin_mem->kin_globalstrategy = globalstrategy;
  kin_mem->kin_constraints = constraints;
  kin_mem->kin_uscale = uscale;
  kin_mem->kin_fscale = fscale;
  kin_mem->kin_precondflag = FALSE;  /* set to the correct state in KINSpgmr */

  if (printfl>0) fprintf(errfp,"scsteptol used: %12.3g \n",scsteptol);
  if (printfl>0) fprintf(errfp,"fnormtol  used: %12.3g \n",fnormtol);

  /*  Calculate the default value for mxnewtstep (max Newton step). */
  if (mxnewtstep == ZERO)
    mxnewtstep = THOUSAND * N_VWL2Norm(uu,uscale);
  if (mxnewtstep < ONE) mxnewtstep = ONE;

  /* Set up the coefficients for the eta calculation. */
  callForcingTerm = ( etaflag != ETACONSTANT);

  /* This value ALWAYS used for Choice 1. */
  if (etaflag == ETACHOICE1)
    ealpha = ONE + HALF * RSqrt(FIVE);

  /* Initial value for eta set to 0.5 for other than the ETACONSTANT option. */
  if (etaflag != ETACONSTANT) 
    eta = HALF;

  /* Initialize all the counters. */
  kin_mem->kin_nfe     = 0;
  kin_mem->kin_nnilpre = 0;
  kin_mem->kin_nni     = 0;
  kin_mem->kin_nbcf    = 0;
  kin_mem->kin_nbktrk  = 0;

  /* See if the system func(uu) = 0 is satisfied by the initial guess uu. */
  if (KINInitialStop(kin_mem)) return(KINSOL_INITIAL_GUESS_OK);

  /* Initialize the linear solver */
  ret = linit(kin_mem);
  if (ret!=0) {
    fprintf(errfp, MSG_LINIT_FAIL);
    return(KINSOL_INPUT_ERROR);
  }

  /* Initialize the L2 norms of f for the linear iteration steps.  */
  fnorm = N_VWL2Norm(fval, fscale);
  f1norm = HALF * fnorm * fnorm;

  if (printfl > 0) 
    fprintf(errfp, "KINSolInit nni= %4d  fnorm= %26.16g  nfe=%6d \n",
            nni, fnorm, nfe);

  /* Problem has been successfully initialized */

  return(SUCCESS);
}

/***************** KINConstraint*************************************
  
 This routine checks the proposed new step uu + x.  If a constraint
 has been violated, then it calculates a new size stepl to be used in
 the global strategy routines (KINInexactNewton or KINLineSearch).

********************************************************************/

static int KINConstraint(KINMem kin_mem) 
{
 realtype mxchange;

 N_VLinearSum(ONE,uu,ONE,pp,vtemp1);
 
 /* The next NVECTOR routine call returns TRUE if all products v1[i]*v2[i]
    are positive (with the proviso that all products which would result
    from v1[i] = 0 are ignored), and FALSE otherwise (e.g. at least one
    such product is negative). */

 if(N_VConstrProdPos(constraints,vtemp1)) return(0);
 
 N_VDiv(pp,uu,vtemp2);
 mxchange = N_VMaxNorm(vtemp2);
 
 if(mxchange>= relu){
   stepl = POINT9 * relu / mxchange;
   return(1);
 }
 return(0);
}


/***************** KINForcingTerm ******************************************
  
 This routine computes eta, the scaling factor in the linear convergence
 stopping tolerance eps when Choice 1 or Choice 2 forcing terms are used.
 Eta is computed here for all but the first iterative step, which is set
 to the default in routine KINSolInit.

 This routine written by Homer Walker of Utah State Univ. with subsequent 
 modifications by Allan G. Taylor @LLNL.

 It is based on the concepts of the paper 'Choosing the forcing terms in 
 an inexact Newton method', SIAM J Sci Comput, 17 (1996), pp 16 - 32,  
 or Utah State University Research Report 6/94/75 of the same title.

 The file kinsol.h documents the input parameters ealpha and egamma as well
 as the selection input parameter for the method of computing eta to be used 
 (that is, Choice 1, Choice 2, and constant eta).  See notes with the inputs 
 iopt[ETACHOICE], ropt[ETAGAMMA], ropt[ETAALPHA], and ropt[ETACONST] 


****************************************************************************/

static void KINForcingTerm(KINMem kin_mem, realtype fnormp)
{
  realtype eta_max = POINT9, eta_min = POINTOHOHOHONE, 
           eta_safe=0.5, linmodel_norm; 

  if (etaflag == ETACHOICE1) {      /* Choice 1 forcing terms. */

    /* Compute the norm of f + J p , scaled L2 norm.  */

    linmodel_norm = RSqrt(fnorm*fnorm + TWO*sfdotJp + sJpnorm*sJpnorm);

    /* Form the safeguarded Choice 1 eta. */ 

    eta_safe = RPowerR(eta, ealpha); 
    eta = ABS(fnormp - linmodel_norm)/fnorm; 

  }

  if (etaflag == ETACHOICE2) {      /* Choice 2 forcing terms. */

    eta_safe = egamma*RPowerR(eta, ealpha); 
    eta = egamma*RPowerR(fnormp/fnorm, ealpha); 

  }

  /* Apply the safeguards. */ 
  if(eta_safe < POINT1) eta_safe = ZERO;
  eta = MAX(eta, eta_safe); 
  eta = MAX(eta, eta_min); 
  eta = MIN(eta, eta_max); 

  return; 

}


/****************** KINFreeVectors *********************************
  
 This routine frees the KINSol vectors allocated in KINAllocVectors.

********************************************************************/

static void KINFreeVectors(KINMem kin_mem)
{
  N_VFree(unew );
  N_VFree(fval );
  N_VFree(pp );
  N_VFree(vtemp1 );
  N_VFree(vtemp2 );
}


/*************** KINInitialStop **********************************

 This routine checks the initial guess (uu) to see if the
 system  func(uu) = 0 is satisfied to the level of 0.01 * fnormtol

******************************************************************/

static booleantype KINInitialStop(KINMem kin_mem)
{
  realtype fmax;

  func(uu, fval, f_data);    nfe++;
  fmax = KINScFNorm(fval, fscale, vtemp1);
  if (printfl > 1) fprintf(errfp,
                 "KINInitialStop: scaled f norm (for stopping)= %12.3g\n",fmax);
  return(fmax <= POINTOH1*fnormtol);
}


/***************************KINInexactNewton**************************

 This routine is the main driver for the inexact Newton algorithm.
 Its purpose is to compute unew = uu + pp in the direction pp from uu,
 taking the full inexact Newton  step.  The step may be constrained if
 the constraint conditions are violated, or if the norm of pp is
 greater than mxnewtstep. 

**********************************************************************/

static int KINInexactNewton(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                            booleantype *maxStepTaken)
{
  int ret;
  realtype pnorm, ratio, ratio1;

  *maxStepTaken = FALSE;
  pnorm = N_VWL2Norm(pp,uscale);
  ratio = ONE;
  if (pnorm > mxnewtstep) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio,pp,pp);
    pnorm = mxnewtstep;
  }

  stepl = pnorm;
  if (printfl > 0) fprintf(errfp,
         " ----- in routine KINInexactNewton (pnorm= %12.4e ) -----\n",pnorm);

  /* If constraints are active, constrain the step accordingly. */
       
  if (constraintsSet) {
  loop{
    ret = KINConstraint(kin_mem);  /* NOTE: This routine changes the step pp. */
    if (ret == 1) {
      ratio1 = stepl / pnorm;
      ratio *= ratio1;
      N_VScale(ratio1, pp, pp);
      pnorm = stepl;
      if (printfl > 0) fprintf(errfp,
		   " --- in routine KINInexactNewton (pnorm= %12.4e \n", pnorm);
      if (pnorm <= scsteptol) return(1);
    }
    else
      break;
    }
  }
 
  /* Scale these two expressions by ratio for later use in KINForcingTerm. */
  sfdotJp = sfdotJp * ratio;
  sJpnorm = sJpnorm * ratio;
 
  /* Compute the iterate unew = uu + pp. */
  N_VLinearSum(ONE, uu, ONE, pp, unew);

  /* Evaluate func(unew) and its norm, and return. */ 
  func(unew, fval, f_data);   nfe++;
  *fnormp = N_VWL2Norm(fval,fscale);
  *f1normp = HALF * (*fnormp) * (*fnormp);
  if (printfl > 1) fprintf(errfp," fnorm (L2) = %20.8e\n", (*fnormp));
  if (pnorm > POINT99 * mxnewtstep ) *maxStepTaken = TRUE; 
  return(0);

}


/***************************************************************************
 The routine KINLineSearch implements the LineSearch algorithm.  Its purpose
 is to find a new unew = uu + rl * pp in the direction pp from uu so that
                                      t
    func(unew) <= func(uu) + alpha * g  (unew - uu)   (alpha = 1.e-4)

        and
                                     t
    func(unew) >= func(uu) + beta * g  (unew - uu)   (beta = 0.9)

    where 0 < rl <= 1. 
****************************************************************************/

static int KINLineSearch(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  int ret, ivio, nfesav, rladjust=0;
  realtype pnorm, ratio, ratio1, slpi, rlmin, rlength, rl, rlmax, rldiff;
  realtype rltmp, rlprev,pt1trl, f1nprv, f1lo, rllo, rlincr, alpha, beta;
  realtype alpha_cond, beta_cond;

  *maxStepTaken = FALSE;
  pnorm = N_VWL2Norm(pp,uscale);
  ratio = ONE; alpha = POINTOHOHOHONE; beta = POINT9;
  if(pnorm > mxnewtstep ){
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio,pp,pp);
    pnorm = mxnewtstep;
  }

  stepl = pnorm;
  ivio = 0;

  /* Check if constraints are active; constrain the step by the constraints. */
  if(constraintsSet){
  loop{
    ret=KINConstraint(kin_mem);
    if(ret == 1){
      ivio = 1;
      ratio1 = stepl / pnorm;
      ratio *= ratio1;
      N_VScale(ratio1,pp,pp);
      pnorm = stepl;

      if(printfl > 0)fprintf(errfp," --- in routine KINLineSearch"
        "(ivio=1, pnorm= %12.4e )\n",pnorm);
    }
    else
     break;
    }
  }

  slpi = sfdotJp * ratio;
  rlength=KINScSteplength(kin_mem, uu, pp, uscale);
  rlmin = scsteptol / rlength;
  rl = ONE;

  if (printfl > 2) fprintf(errfp,"KINLineSearch -----\n"
       "  min_lam=%11.4e  f1norm=%11.4e  pnorm=%11.4e\n",rlmin, f1norm, pnorm);

  /* Now begin the iteration to find a rl value which satisfies both the
     alpha- and beta- conditions. If rl < rlmin, terminate and return 1.  */

  nfesav = nfe;

  loop {

    N_VLinearSum(ONE, uu, rl, pp, unew);
    func(unew, fval, f_data);   nfe++;
    *fnormp = N_VWL2Norm(fval,fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp) ;
    alpha_cond = f1norm + alpha*slpi*rl;
    if (printfl > 2 && rladjust > 0) fprintf(errfp,
         " fnormp=%15.8e  f1normp=%15.8e  alpha_cond=%15.8e lam=%15.8e\n", 
         *fnormp,*f1normp, alpha_cond, rl);
    if ((*f1normp) <= alpha_cond) break;

    /* alpha condition not satisfied; perform quadratic backtrack to compute
       new rl value.  */
    if (rl < rlmin) {
      /* No satisfactory unew can be found sufficiently distinct from uu.
         Copy uu into unew and return.  */
      N_VScale(ONE, uu, unew);
      func(unew, fval, f_data);   nfe++;
      *fnormp = N_VWL2Norm(fval,fscale);
      *f1normp = HALF * (*fnormp) * (*fnormp);
      nbktrk = nfe - nfesav - 1;
      return(1);
    }
    rltmp = -slpi / (TWO*((*f1normp)-f1norm-slpi));
    if (rltmp > HALF * rl) rltmp = HALF * rl;
    rlprev = rl;
    f1nprv = (*f1normp);
    pt1trl = POINT1*rl;
    rl = MAX(pt1trl,rltmp);
    rladjust++;
  }

  /* The alpha condition is satisfied. Now check the beta condition. */

  beta_cond = f1norm + beta*slpi*rl;
  if (*f1normp < beta_cond) {
    if (rl == ONE && pnorm < stepl) {
      rlmax = stepl / pnorm;
      if (ivio == 1) rlmax = ONE;
      do {
        rlprev = rl;
        f1nprv = *f1normp;
        rl = MIN(TWO*rl,rlmax);
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data);   nfe++;
        *fnormp = N_VWL2Norm(fval,fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);
        alpha_cond = f1norm + alpha*slpi*rl;
        beta_cond = f1norm + beta*slpi*rl;
        if (printfl > 2) fprintf(errfp,"  f1normp=%15.8e  beta_cond=%15.8e  "
                         "lam=%15.8e\n", *f1normp, beta_cond, rl);
      } while ((*f1normp) <= alpha_cond && 
          (*f1normp) < beta_cond && (rl < rlmax) );
    }
    if (rl < ONE || ((rl > ONE) && (*f1normp > alpha_cond) ) ) {
      rllo = MIN(rl,rlprev);
      rldiff = ABS(rlprev - rl);

      if (rl < rlprev) f1lo = *f1normp;
      else f1lo = f1nprv;

      do {
        rlincr = HALF * rldiff;
        rl = rllo + rlincr;
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data);   nfe++;
        *fnormp = N_VWL2Norm(fval,fscale);
        *f1normp = HALF * *fnormp * *fnormp;
        alpha_cond = f1norm + alpha*slpi*rl;
        beta_cond = f1norm + beta*slpi*rl;
        if (printfl > 2&& rladjust>0)
        fprintf(errfp,"  f1normp=%12.5e  alpha_cond=%12.5e  beta_cond=%12.5e"
                "  lam=%12.5e\n", *f1normp, alpha_cond, beta_cond, rl);

        if ((*f1normp) > alpha_cond) 
          rldiff = rlincr;
        else if (*f1normp < beta_cond) {
          rllo = rl;
          rldiff = rldiff - rlincr;
          f1lo = *f1normp;
        }

      } while (*f1normp > alpha_cond ||
        ((*f1normp < beta_cond) && (rldiff > rlmin)));

      if ((*f1normp) < beta_cond) {
      /*  beta condition could not be satisfied; set unew to last u value
          that satisfied the alpha condition and continue.  Increment
          counter on number of steps beta condition not satisfied.       */

        N_VLinearSum(ONE, uu, rllo, pp, unew);
        func(unew, fval, f_data);   nfe++;
        *fnormp = N_VWL2Norm(fval,fscale);
        *f1normp = HALF * *fnormp * *fnormp;   
        nbcf++;
      }
    }  /* end of if (rl < ONE) block */
  } /* end of (f1normp < test) loop */

  nbktrk= nfe - nfesav - 1;
  if (printfl > 1 && rladjust>0) fprintf(errfp, "Number of lambda adjustments ="
            "%d\n", rladjust);

  /* Scale these two expressions by rl*ratio for later use in KINForcingTerm. */
  sfdotJp = sfdotJp * rl * ratio;
  sJpnorm = sJpnorm * rl * ratio;

  if (rl*pnorm > POINT99 * mxnewtstep) *maxStepTaken = TRUE;
  return(0);
}
  
/**************************KINLinSolDrv*************************************

 This routine handles the process of solving for the approximate solution
 xx of the Newton equations in the Newton iteration.  Subsequent routines 
 handle the nonlinear aspects of its applications. It takes as arguments 
 the KINSol memory pointer, an N_Vector, called xx, which is used as scratch 
 and when returned contains the solution x!  It also takes an N_Vector 
 bb which is used as scratch but represents the rhs of the system J x = b
 whose solution xx is being approximated.

****************************************************************************/

static int KINLinSolDrv(KINMem kin_mem , N_Vector bb , N_Vector xx )
{
  int ret;

  if (nni-nnilpre >= msbpre) pthrsh = TWO;

  loop{

    precondcurrent = FALSE;

    if( (pthrsh > ONEPT5) && setupNonNull ) {
      ret = lsetup(kin_mem);  /* Call precondset routine.  */
      precondcurrent = TRUE;
      nnilpre = nni;
      if(ret != 0) return(KINSOL_PRECONDSET_FAILURE);
    }

  /*  Load bb with the current value of -fval.  */
             
  N_VScale(-ONE, fval, bb);

  /*  Call the generic 'lsolve' routine to solve the system J x = b.  */

  ret = lsolve(kin_mem, xx, bb, &res_norm);

  if (ret != 1)return(ret);

  if (!precondflag)return(KINSOL_KRYLOV_FAILURE);

  if (precondcurrent)return(KINSOL_KRYLOV_FAILURE);

  /* Loop back only if the preconditioner is in use and is not current. */

  pthrsh = TWO;

  }  /*  end of loop  */

}


/***********************KINScFNorm *********************************

 This routine computes the max norm for scaled vectors. The scaling
 vector is scale, the vector of which the norm is to be determined
 is vv. The returned value, fnormval, is the resulting scaled vector
 norm.  This routine uses N_Vector functions from the vector module.

********************************************************************/

static realtype KINScFNorm(N_Vector vv , N_Vector scale , N_Vector wrkv)
{
  N_VAbs(vv, wrkv);
  N_VProd(scale, wrkv, wrkv);
  return(N_VMaxNorm(wrkv));
}

/************************KINScSteplength ***************************

 This routine computes the max norm of the scaled steplength, ss.
 Here ucur is the current step and usc is the u scale factor.

********************************************************************/

static realtype KINScSteplength(KINMem kin_mem, N_Vector ucur,
                                N_Vector ss, N_Vector usc)
{
  N_VInv(usc, vtemp1);
  N_VAbs(ucur, vtemp2);
  N_VLinearSum(ONE, vtemp1, ONE, vtemp2, vtemp1);
  N_VDiv(ss, vtemp1, vtemp1);
  return(N_VMaxNorm(vtemp1));
}

/************** KINStop *****************************************

 This routine checks the current iterate unew to see if the
 system  func(unew) = 0 is satisfied by a variety of tests.

*****************************************************************/

static int KINStop(KINMem kinmem, booleantype maxStepTaken, int globalstratret)
{
  KINMem kin_mem;
  realtype fmax, rlength;

  kin_mem=kinmem;
  if (globalstratret==1){
    if (precondflag && setupNonNull && !precondcurrent) {
      /*  If the globalstratret was caused (potentially) by the 
          preconditioner being out of date, update the preconditioner. */
      pthrsh = TWO;
      return(CONTINUE_ITERATIONS);
    }
    else return( (globalstrategy == INEXACT_NEWTON) ?
                  KINSOL_STEP_LT_STPTOL : KINSOL_LNSRCH_NONCONV );
  }

  /*  Check tolerance on scaled norm of func at the current iterate. */

  fmax = KINScFNorm(fval, fscale, vtemp1);
  if (printfl>1) fprintf(errfp," scaled f norm (for stopping)= %12.3g\n",fmax);

  if (fmax <= fnormtol) return(SUCCESS);

  /*  Check for the scaled distance between the last two steps too small. */

  N_VLinearSum(ONE,unew,-ONE,uu,vtemp1);
  rlength = KINScSteplength(kin_mem, unew, vtemp1, uscale);
  if (rlength <= scsteptol) {
    if (!precondcurrent) {
      /* For rlength too small and the preconditioner not current,
         try again with the preconditioner current.                */
      pthrsh = TWO;
      return(CONTINUE_ITERATIONS);
    }
    else
      /* Otherwise return failure flag. */
      return(KINSOL_STEP_LT_STPTOL);
  }

 /* If the max. number of iterations is reached, return failure flag. */
   if (nni >= mxiter) return(KINSOL_MAXITER_REACHED);

 /* Check for consecutive number of steps taken of size mxnewtstep.
    If not maxStepTaken, set ncscmx to 0 */
 
  if (maxStepTaken) ncscmx++;
  else ncscmx = 0;
 
  if (ncscmx == 5) return(KINSOL_MXNEWT_5X_EXCEEDED);

  /* Load threshold for re-evaluating the preconditioner. */
  pthrsh = rlength;

  /* If make it to here, the iteration process is not finished; 
     return CONTINUE_ITERATIONS. */
  return(CONTINUE_ITERATIONS);
}


/*******************************************************************/
/********* END Private Helper Functions Implementation *************/
/*******************************************************************/


/****************************************************************/
/************** END KINSol Implementation ***********************/
/****************************************************************/
