/******************************************************************
 *                                                                *
 * File          : senspvode.c                                    *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Last Modified : 05 February 2001                               *
 *----------------------------------------------------------------*
 * This is the implementation file for computing the sensitivity  *
 * of the ODE solution y(t, p) with respect to the parameters     *
 * in p.                                                          *
 *                                                                *
 ******************************************************************/


/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cvode.h"    /* CVODE constants and prototypes                    */
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "llnlmath.h" /* header file for C math library                    */
#include "senspvode.h" /* sensitivity data types and prototypes            */

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/

/************************************************************/
/************** BEGIN Private Constants *********************/
/************************************************************/

#define ZERO      RCONST(0.0)  /* real 0.0   */
#define HALF      RCONST(0.5)  /* real 0.5   */
#define ONE       RCONST(1.0)  /* real 1.0   */

/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* SensCVodeMalloc/SensCVReInit Error Messages */

#define SensCVM         "SensCVodeMalloc/SensCVReInit-- "

#define MSG_REI_NO_MEM   "SensCVReInit-- cvode_mem = NULL illegal.\n\n"
#define MSG_REI_NO_SDATA "SensCVReInit-- s_data = NULL illegal.\n\n"
#define MSG_BAD_NY       SensCVM "Ny=%ld < 1 illegal.\n\n"
#define MSG_BAD_NS       SensCVM "Ns=%ld < 0 illegal.\n\n"
#define MSG_P_NULL       SensCVM "p=NULL illegal.\n\n"
#define MSG_PBAR_NULL    SensCVM "pbar=NULL illegal.\n\n"
#define MSG_F_NULL       SensCVM "f=NULL illegal.\n\n"
#define MSG_FTEMP_NULL   SensCVM "ftemp=NULL illegal.\n\n"
#define MSG_YTEMP_NULL   SensCVM "ytemp=NULL illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/

/************************************************************/
/************** END Private Constants ***********************/
/************************************************************/


/**************************************************************/
/**************** BEGIN Readability Constants *****************/
/**************************************************************/

#define CENTERED1  0
#define CENTERED2  1
#define FORWARD1   2
#define FORWARD2   3
#define nfe        (cvmem->cv_nfe)

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/

/**************************************************************/
/********* BEGIN Helper Function Prototypes *******************/
/**************************************************************/

void FDQ(integer Ntotal, real t, N_Vector u, N_Vector udot, void *s_data);

/**************************************************************/
/********** END Helper Function Prototypes ********************/
/**************************************************************/


/***************************************************************/
/************* BEGIN SENSITIVITY Implementation ****************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/


/******************** SensCVodeMalloc ***************************

 SensCVodeMalloc initializes a block of sensitivity data.
 It also allocates and initializes memory for the ODE problem. 
 All problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.
 
*****************************************************************/

void *SensCVodeMalloc(int Ny, int Ns, int Ntotal, RhsFn f, real t0, 
		      N_Vector y0, int lmm, int iter, int itol, 
		      real *reltol, void *abstol, void *f_data, 
		      FILE *errfp, boole optIn, long int iopt[], 
		      real ropt[], void *machEnv, real p[], 
		      real pbar[], void *plist, real rhomax)
{
  SensData sdata;
  void *cv_mem;
  FILE *fp;

  /* Check for legal input parameters */

  fp = (errfp == NULL) ? stdout : errfp;

  if (Ny < 1) {
    fprintf(fp, MSG_BAD_NY, Ny);
    return(NULL);
  }
  if (Ns < 0) {
    fprintf(fp, MSG_BAD_NS, Ns);
    return(NULL);
  }
  if (f == NULL) {
    fprintf(fp, MSG_F_NULL);
    return(NULL);
  }
  if (p == NULL) {
    fprintf(fp, MSG_P_NULL);
    return(NULL);
  }
  if (pbar == NULL) {
    fprintf(fp, MSG_PBAR_NULL);
    return(NULL);
  }

  sdata = (SensData) malloc(sizeof(*sdata));
  sdata->Ny = Ny;
  sdata->Ns = Ns;
  sdata->f_data = f_data;
  sdata->p = p;
  sdata->pbar = pbar;
  sdata->plist = plist;
  sdata->rhomax = rhomax;
  sdata->f_ptr = f; 
  sdata->ftemp = N_VNew(Ny, machEnv);
  sdata->ytemp = N_VNew(Ny, machEnv);
  cv_mem = CVodeMalloc(Ntotal, FDQ, t0, y0, lmm, iter, itol, reltol,
		       abstol, sdata, errfp, optIn, iopt, ropt, machEnv);
  sdata->cv_mem = cv_mem;
  return(cv_mem);
}


/******************** SensCVReInit **********************************

 SensCVReInit reinitializes CVODE's memory for a problem, assuming
 it has already been allocated in a prior SensCVodeMalloc call.
 All problem specification inputs are checked for errors.
 The problem size parameters (Ny, Ns, Ntotal) are assumed to be 
 unchanged since the call to SensCVodeMalloc, and the maximum order 
 maxord must not be larger.
 If any error occurs during initialization, it is reported to the
 file whose file pointer is errfp.
 The return value is SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
 
*****************************************************************/

int SensCVReInit(void *cvode_mem, RhsFn f, real t0, N_Vector y0,
             int lmm, int iter, int itol, real *reltol, void *abstol,
             void *f_data, FILE *errfp, boole optIn, long int iopt[],
             real ropt[], void *machEnv, real p[], real pbar[],
             void *plist, real rhomax)
{
  CVodeMem cv_mem;
  SensData sdata;
  FILE *fp;
  int flag;

  /* Check for legal input parameters */

  fp = (errfp == NULL) ? stdout : errfp;

  if (cvode_mem == NULL) {
    fprintf(fp, MSG_REI_NO_MEM);
    return(SensCVREI_NO_MEM);
  }

  cv_mem = (CVodeMem) cvode_mem;

  if (p == NULL) {
    fprintf(fp, MSG_P_NULL);
    return(SensCVREI_ILL_INPUT);
  }
  if (pbar == NULL) {
    fprintf(fp, MSG_PBAR_NULL);
    return(SensCVREI_ILL_INPUT);
  }
  if (f == NULL) {
    fprintf(fp, MSG_F_NULL);
    return(NULL);
  }
  if (cv_mem->cv_f_data == NULL) {
    fprintf(fp, MSG_REI_NO_SDATA);
    return(SensCVREI_ILL_INPUT);
  }

  sdata = (SensData) cv_mem->cv_f_data;

  if (sdata->Ny < 1) {
    fprintf(fp, MSG_BAD_NY, sdata->Ny);
    return(SensCVREI_ILL_INPUT);
  }
  if (sdata->Ns < 0) {
    fprintf(fp, MSG_BAD_NS, sdata->Ns);
    return(SensCVREI_ILL_INPUT);
  }
  if (sdata->ftemp == NULL) {
    fprintf(fp, MSG_FTEMP_NULL);
    return(SensCVREI_ILL_INPUT);
  }
  if (sdata->ytemp == NULL) {
    fprintf(fp, MSG_YTEMP_NULL);
    return(SensCVREI_ILL_INPUT);
  }

  sdata->f_ptr = f;
  sdata->f_data = f_data;
  sdata->p = p;
  sdata->pbar = pbar;
  sdata->plist = plist;
  sdata->rhomax = rhomax;

  flag = CVReInit(cv_mem, FDQ, t0, y0, lmm, iter, itol, reltol,
		  abstol, sdata, errfp, optIn, iopt, ropt, machEnv);

  return(flag);
}


/********************* SensCVodeFree ******************************

 This routine frees the block of sensitivity data and problem memory 
 allocated by SensCVodeMalloc. Such memory includes all the vectors 
 allocated by CVAllocVectors, and the memory lmem for the linear 
 solver. 

*******************************************************************/

void SensCVodeFree(void *cvode_mem)
{
  CVodeMem cv_mem;
  SensData sdata;

  cv_mem = (CVodeMem) cvode_mem;
  sdata = (SensData) cv_mem->cv_f_data;

  N_VFree(sdata->ftemp);
  N_VFree(sdata->ytemp);
  free(sdata);

  CVodeFree(cvode_mem);
}


/********************* SensSetVecAtol *****************************

 This routine sets the vector absolute error tolerances for the Ns
 scaled sensitivity vectors to be the same as the vector absolute
 error tolerances for the ODEs contained in y' = f(t,y,p).

*******************************************************************/

void SensSetVecAtol(N_Vector atol, integer Ns)
{
  N_Vector *atolsub;
  int i;

  atolsub = N_VSUB(atol);

  for (i = 1; i <= Ns; i++)
    N_VScale(ONE, atolsub[0], atolsub[i]);
}


/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Helper Function Implementation *********************/
/*******************************************************************/

/****************** FDQ ********************************************

 This routine evaluates the right hand side of the ODE 

      y' = f(t,y,p)

 and approximates the right hand side of the scaled sensitivity ODE

      w'_i = (df/dy)*w_i + pbar_i*(df/dp_i)

 using finite differences.

**********************************************************************/

void FDQ(integer Ntotal, real t, N_Vector u, N_Vector udot, void *s_data)
{
  integer i, j, Ns, Ny, method;
  int *plist;
  real pbar, psave;
  real delta, rdelta, r2delta;
  real deltap, deltay, rdeltap, rdeltay, r2deltap, r2deltay;
  real rtol, uround, normw;
  real ratio, rhomax;
  N_Vector *usub, *udotsub, *ewtsub;
  N_Vector ftemp, ytemp, vtemp;
  CVodeMem cvmem;
  SensData sdata;
  void *fdata;

  sdata = (SensData) s_data;
  cvmem = (CVodeMem) sdata->cv_mem;

  usub = N_VSUB(u);
  udotsub = N_VSUB(udot);
  ewtsub = N_VSUB(cvmem->cv_ewt);

  /* Extract parameters and workspace from sdata */
  Ns = sdata->Ns;
  Ny = sdata->Ny;
  plist = sdata->plist;
  fdata = sdata->f_data;
  ftemp = sdata->ftemp;
  ytemp = sdata->ytemp;

  rtol = *(cvmem->cv_reltol);
  uround = cvmem->cv_uround;

  /******************************************/
  /* Set deltap for perturbing a parameter  */
  /******************************************/
  deltap = RSqrt(MAX(rtol, uround));
  rdeltap = ONE/deltap;

  /*******************************************************************
   * Compute y derivative as y' = f(t,y,p).                          *
   *******************************************************************/
  sdata->f_ptr(Ny, t, usub[0], udotsub[0], fdata);
  
  for (i = 1; i <= Ns; i++) {

    j = (plist == NULL) ? i : plist[i-1];

  /*******************************************************************
   * Estimate i-th scaled sensitivity derivative w'_i using          *
   * finite difference formulas.                                     *
   *******************************************************************/

    psave = sdata->p[j-1];
    pbar = sdata->pbar[j-1];

    /*************************************/
    /* Determine deltay for perturbing y */
    /*************************************/
    normw = N_VWrmsNorm(usub[i], ewtsub[0]);
    rdeltay = MAX(normw, rdeltap);
    deltay = ONE/rdeltay;

    /***********************************************************/
    /* Select finite difference formula for w'_i               */
    /***********************************************************/
    ratio = deltay*rdeltap;
    rhomax = sdata->rhomax;

    if ((MAX(ONE/ratio, ratio) <= ABS(rhomax)) || rhomax == ZERO)
      /***************************************************/
      /* Centered or forward difference formula for w'_i */
      /***************************************************/
      method = (rhomax >= ZERO) ? CENTERED1 : FORWARD1; 
    else {
      /*****************************************************************/
      /* Centered or forward differences for each term of w'_i:        */
      /* (df/dy)*w_i and pbar_i*(df/dp_i).                             */
      /*****************************************************************/
      method = (rhomax > ZERO) ? CENTERED2 : FORWARD2;
    }

    switch (method) {

      /**********************************************/
    case 0: /* Centered difference formula for w'_i */
      /**********************************************/
    
      delta = MIN(deltay, deltap);
      r2delta = HALF/delta;

      /* Forward perturb y and parameter */
      N_VLinearSum(delta, usub[i], ONE, usub[0], ytemp);
      sdata->p[j-1] = psave + delta*pbar;

      sdata->f_ptr(Ny, t, ytemp, udotsub[i], fdata);
      nfe++;

      /* Backward perturb y and parameter */
      N_VLinearSum(-delta, usub[i], ONE, usub[0], ytemp);
      sdata->p[j-1] = psave - delta*pbar;

      sdata->f_ptr(Ny, t, ytemp, ftemp, fdata);
      nfe++;

      /* Centered differencing for 2nd order accuracy */
      N_VLinearSum(r2delta, udotsub[i], -r2delta, ftemp, udotsub[i]);
  
      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

      /**************************************************************/
    case 1: /* Centered differences for each term of w'_i:          */
            /* (df/dy)*w_i and pbar_i*(df/dp_i)                     */       
      /**************************************************************/

      r2deltap = HALF/deltap;
      r2deltay = HALF/deltay;
      
      /* Forward perturb y only */
      N_VLinearSum(deltay, usub[i], ONE, usub[0], ytemp);

      sdata->f_ptr(Ny, t, ytemp, udotsub[i], fdata);
      nfe++;

      /* Backward perturb y only */
      N_VLinearSum(-deltay, usub[i], ONE, usub[0], ytemp);

      sdata->f_ptr(Ny, t, ytemp, ftemp, fdata);
      nfe++;

      /*************************************/
      /* Store (df/dy)*w_i into udotsub[i] */
      /*************************************/
      N_VLinearSum(r2deltay, udotsub[i], -r2deltay, ftemp, udotsub[i]);

      /* Rename ytemp as vtemp for readability */
      vtemp = ytemp;

      /* Forward perturb parameter */
      sdata->p[j-1] = psave + deltap*pbar;
      
      sdata->f_ptr(Ny, t, usub[0], vtemp, fdata);
      nfe++;
      
      /* Backward perturb parameter */
      sdata->p[j-1] = psave - deltap*pbar;

      sdata->f_ptr(Ny, t, usub[0], ftemp, fdata);
      nfe++;

      /***************************************/
      /* Store pbar_i*(df/dp_i) into vtemp   */
      /***************************************/
      N_VLinearSum(r2deltap, vtemp, -r2deltap, ftemp, vtemp);

      /*** Compute udotsub[i] = (df/dy)*w_i + pbar_i*(df/dp_i) ***/
      N_VLinearSum(ONE, udotsub[i], ONE, vtemp, udotsub[i]);

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

      /**********************************************/
    case 2: /* Forward difference formula for w'_i  */
      /**********************************************/
    
      delta = MIN(deltay, deltap);
      rdelta = ONE/delta;

      /* Forward perturb y and parameter */
      N_VLinearSum(delta, usub[i], ONE, usub[0], ytemp);
      sdata->p[j-1] = psave + delta*pbar;

      sdata->f_ptr(Ny, t, ytemp, ftemp, fdata);
      nfe++;

      /* Forward differencing for 1st order accuracy */
      N_VLinearSum(rdelta, ftemp, -rdelta, udotsub[0], udotsub[i]);
  
      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

      /*************************************************************/
    case 3: /* Forward differences for each term of w'_i:          */
            /* (df/dy)*w_i and pbar_i*(df/dp_i).                   */
      /*************************************************************/

      /* Perturb y only*/
      N_VLinearSum(deltay, usub[i], ONE, usub[0], ytemp);

      sdata->f_ptr(Ny, t, ytemp, ftemp, fdata);
      nfe++;

      /*************************************/
      /* Store (df/dy)*w_i into udotsub[i] */
      /*************************************/
      N_VLinearSum(rdeltay, ftemp, -rdeltay, udotsub[0], udotsub[i]);

      /* Forward perturb the parameter */
      sdata->p[j-1] = psave + deltap*pbar;

      sdata->f_ptr(Ny, t, usub[0], ftemp, fdata);
      nfe++;

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      /* Rename ytemp as vtemp for readability */
      vtemp = ytemp;

      /***************************************/
      /* Store pbar_i*(df/dp_i) into vtemp   */
      /***************************************/
      N_VLinearSum(rdeltap, ftemp, -rdeltap, udotsub[0], vtemp);

      /*** Compute udotsub[i] = (df/dy)*w_i + pbar_i*(df/dp_i) ***/
      N_VLinearSum(ONE, udotsub[i], ONE, vtemp, udotsub[i]);

      break;

    default:      
      /* do nothing */
      break;
    }
  }
}


/*******************************************************************/
/******** END Helper Function Implementation ***********************/
/*******************************************************************/


/***************************************************************/
/************** END SENSITIVITY Implementation *****************/
/***************************************************************/
