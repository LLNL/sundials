/******************************************************************
 *                                                                *
 * File          : sensida.c                                      *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 05 February 2001                               *
 *----------------------------------------------------------------*
 * This is the implementation file for computing the sensitivity  *
 * of the DAE solution y(t,p) with respect to the parameters      *
 * in p.                                                          *
 *                                                                *
 ******************************************************************/


/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ida.h"      /* IDA constants and prototypes                      */
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "llnlmath.h" /* header file for C math library                    */
#include "sensida.h"  /* sensitivity data types and prototypes             */

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/

/************************************************************/
/************** BEGIN Private Constants *********************/
/************************************************************/

#define ZERO      RCONST(0.0)  /* real 0.0   */
#define HALF      RCONST(0.5)  /* real 0.5   */
#define ONE       RCONST(1.0)  /* real 1.0   */
#define TWO       RCONST(2.0)  /* real 2.0   */

/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* SensIDAMalloc error messages */

#define SensIDAM           "SensIDAMalloc-- "

#define MSG_BAD_NY         SensIDAM "Ny=%ld < 1 illegal.\n\n"
#define MSG_BAD_NS         SensIDAM "Ns=%ld < 1 illegal.\n\n"
#define MSG_P_NULL         SensIDAM "p=NULL illegal.\n\n"
#define MSG_PBAR_NULL      SensIDAM "pbar=NULL illegal.\n\n"
#define MSG_RES_NULL       SensIDAM "res=NULL illegal.\n\n"

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

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/

/**************************************************************/
/********* BEGIN Helper Function Prototypes *******************/
/**************************************************************/

int RESDQ(integer Ntotal, real t, N_Vector yy, N_Vector yp, 
	  N_Vector resval, void *s_data);

/**************************************************************/
/********** END Helper Function Prototypes ********************/
/**************************************************************/


/***************************************************************/
/************* BEGIN SENSITIVITY Implementation ****************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/


/******************** SensIDAMalloc ***************************

 SensIDAMalloc initializes a block of sensitivity data.
 It also allocates and initializes memory for the DAE problem.
 All problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.

*****************************************************************/

void *SensIDAMalloc(integer Ny, integer Ns, integer Ntotal, ResFn res,
		    void *rdata, real t0, N_Vector y0, N_Vector yp0, 
		    int itol, real *rtol, void *atol, N_Vector id, 
		    N_Vector constraints, FILE *errfp, boole optIn, 
		    long int iopt[], real ropt[], void *machEnv,
		    real p[], real pbar[], void *plist, real rhomax) 
{ 
  void *ida_mem;
  SensData sdata;
  FILE *fp;
  integer retval;

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
  if (p == NULL) {
    fprintf(fp, MSG_P_NULL);
    return(NULL);
  }
  if (pbar == NULL) {
    fprintf(fp, MSG_PBAR_NULL);
    return(NULL);
  }
  if (res == NULL) {
    fprintf(fp, MSG_RES_NULL);
    return(NULL);
  }

  sdata = (SensData) malloc(sizeof(*sdata));
  sdata->Ny = Ny;
  sdata->Ns = Ns;
  sdata->p = p;
  sdata->pbar = pbar;
  sdata->plist = plist;
  sdata->rhomax = rhomax;
  sdata->res_ptr = res;
  sdata->res_data = rdata;
  sdata->restemp = N_VNew(Ny, machEnv);
  sdata->ytemp = N_VNew(Ny, machEnv);
  sdata->yptemp = N_VNew(Ny, machEnv);

  ida_mem = IDAMalloc(Ntotal, RESDQ, sdata, t0, y0, yp0, itol, rtol, atol,
                  id, constraints, NULL, optIn, iopt, ropt, machEnv);
  sdata->ida_mem = ida_mem;

  return(ida_mem);
}

/********************* SensIDAFree ******************************

 This routine frees the block of sensitivity data and problem memory
 allocated by SensIDAMAlloc. Such memory includes all the vectors
 allocated by IDAAllocVectors, and the memory lmem for the linear
 solver.

*******************************************************************/

void SensIDAFree(void *ida_mem)
{
  IDAMem IDA_mem;
  SensData sdata;

  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem == NULL) return;

  sdata = (SensData) IDA_mem->ida_rdata;
  N_VFree(sdata->restemp);
  N_VFree(sdata->ytemp);
  N_VFree(sdata->yptemp);
  free(sdata);

  IDAFree(IDA_mem);
}

/********************* SensSetVecAtol *****************************
 This routine sets the vector absolute error tolerances to the Ns
 scaled sensitivity vectors to be the same as the vector absolute
 error tolerances for the equations contained in F(t,y,y',p) = 0.
*******************************************************************/

void SensSetVecAtol(N_Vector atol, integer Ns)
{
  N_Vector *atolsub;
  int i;

  atolsub = N_VSUB(atol);

  for (i = 1; i <= Ns; i++)
    N_VScale(ONE, atolsub[0], atolsub[i]);
}


/********************* SensSetId **********************************
 SensSetId identifies each component of the Ns sensitivity vectors
 as either a differential or algebraic variable.
 A call to this subroutine is required if the optional input
 SUPPRESSALG is set, or if IDACalcIC is to be called with 
 icopt = CALC_YA_YDP_INIT.
*******************************************************************/

void SensSetId(N_Vector id, integer Ns)
{
  N_Vector *idsub;
  int i;

  idsub = N_VSUB(id);

  for (i = 1; i <= Ns; i++)
    N_VScale(ONE, idsub[0], idsub[i]);

}

/********************* SensInitZero *******************************
 This routine initializes the trailing Ny*Ns components of a vector
 to zero.
*******************************************************************/

void SensInitZero(N_Vector z, integer Ns)
{
  N_Vector *zsub;
  int i;

  zsub = N_VSUB(z);

  for (i = 1; i <= Ns; i++)
    N_VConst(ZERO, zsub[i]);

}


/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Helper Function Implementation *********************/
/*******************************************************************/

/****************** RESDQ ********************************************

 This routine evaluates the residual of the DAE 

     F(t,y,y',p) = 0

 and approximates the residual of the DAE sensitivity equations

     (dF/dy)*w_i + (dF/dy')*w'_i + pbar_i*(dF/dp_i)

 using finite differences.

**********************************************************************/

int RESDQ(integer Ntotal, real t, N_Vector yy, N_Vector yp, 
	  N_Vector resval, void *s_data)
{
  integer i, j, Ny, Ns, method;
  int retval;
  int *plist;
  real pbar, psave;
  real cj;
  real delta, rdelta, r2delta;
  real deltap, rdeltap, r2deltap;
  real deltay, rdeltay, r2deltay;
  real normw;
  real rtol, uround;
  N_Vector *yysub, *ypsub, *resvalsub, *ewtsub;
  N_Vector restemp, ytemp, yptemp;
  IDAMem idamem;
  SensData sdata;
  void *rdata;
  real *resdata;
  real ratio, rhomax;

  sdata = (SensData) s_data;
  idamem = (IDAMem) sdata->ida_mem;

  yysub = N_VSUB(yy);
  ypsub = N_VSUB(yp);
  resvalsub = N_VSUB(resval);
  ewtsub = N_VSUB(idamem->ida_ewt);

  /* extract parameters and workspace from s_data */
  Ns = sdata->Ns;
  Ny = sdata->Ny;
  plist = sdata->plist;
  rdata = sdata->res_data;
  restemp = sdata->restemp;
  ytemp = sdata->ytemp;
  yptemp = sdata->yptemp;

  rtol = *(idamem->ida_rtol);
  uround = idamem->ida_uround;
  cj = idamem->ida_cj;

  /* Set deltap for perturbing a parameter */
  deltap = RSqrt(MAX(rtol, uround));
  rdeltap = ONE/deltap;
  r2deltap = HALF/deltap;

  /* Evaluate the base value of the DAE residual */
  retval = sdata->res_ptr(Ny, t, yysub[0], ypsub[0], resvalsub[0], rdata);
  if (retval != 0) return(retval);
  
  for (i = 1; i <= Ns; i++) {

    j = (plist == NULL) ? i: plist[i-1];

    /* Estimate the residual for the i-th scaled sensitivity equation */

    psave = sdata->p[j-1];
    pbar = sdata->pbar[j-1];

    /* Set deltay for perturbing y and y' */
    normw = N_VWrmsNorm(yysub[i], ewtsub[0]);
    rdeltay = MAX(normw, rdeltap);
    deltay = ONE/rdeltay;
    r2deltay = HALF/deltay;

    /* Select finite difference formula */
    ratio = deltay*rdeltap;
    rhomax = sdata->rhomax;

    if ((MAX(ONE/ratio, ratio) <= ABS(rhomax)) || rhomax == ZERO)
      method = (rhomax >= ZERO) ? CENTERED1 : FORWARD1;
    else {
      method = (rhomax > ZERO) ? CENTERED2 : FORWARD2;
    }

    switch (method) {

      /*********************************************/
    case 0: /* Centered difference formula         */
      /*********************************************/

      delta = MIN(deltay, deltap);
      rdelta = ONE/delta;
      r2delta = HALF/delta;

      /* Forward perturb y, y' and parameter */
      N_VLinearSum(delta, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(delta, ypsub[i], ONE, ypsub[0], yptemp);
      sdata->p[j-1] = psave + delta*pbar;

      /* Save residual in resvalsub[i] */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, resvalsub[i], rdata);
      if (retval != 0) return(retval);

      /* Backward perturb y, y' and parameter */
      N_VLinearSum(-delta, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(-delta, ypsub[i], ONE, ypsub[0], yptemp);
      sdata->p[j-1] = psave - delta*pbar;

      /* Save residual in restemp */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, restemp, rdata);
      if (retval != 0) return(retval);

      /* Estimate the residual for the i-th sensitivity equation */
      N_VLinearSum(r2delta, resvalsub[i], -r2delta, restemp, resvalsub[i]);

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

    break;

    /*************************************************************/
    case 1: /* Centered differences on separate terms:           */
            /* (dF/dy) w_i + (dF/dy')* w'_i + pbar_i * (dF/dp_i) */
      /***********************************************************/
      
      /* Forward perturb y and y' */
      N_VLinearSum(deltay, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(deltay, ypsub[i], ONE, ypsub[0], yptemp);

      /* Save residual in resvalsub[i] */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, resvalsub[i], rdata);
      if (retval != 0) return(retval);

      /* Backward perturb y and y' */
      N_VLinearSum(-deltay, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(-deltay, ypsub[i], ONE, ypsub[0], yptemp);

      /* Save residual in restemp */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, restemp, rdata);
      if (retval != 0) return(retval);

      /* Save the first difference quotient in resvalsub[i] */
      N_VLinearSum(r2deltay, resvalsub[i], -r2deltay, restemp, resvalsub[i]);

      /* Forward perturb parameter */
      sdata->p[j-1] = psave + deltap*pbar;

      /* Save residual in ytemp */
      retval = sdata->res_ptr(Ny, t, yysub[0], ypsub[0], ytemp, rdata);
      if (retval != 0) return(retval);

      /* Backward perturb parameter */
      sdata->p[j-1] = psave - deltap*pbar;

      /* Save residual in yptemp */
      retval = sdata->res_ptr(Ny, t, yysub[0], ypsub[0], yptemp, rdata);
      if (retval != 0) return(retval);

      /* Save the second difference quotient in restemp */
      N_VLinearSum(r2deltap, ytemp, -r2deltap, yptemp, restemp);


      /* Add the difference quotients for the residual for the i-th
         sensitivity */
      N_VLinearSum(ONE, resvalsub[i], ONE, restemp, resvalsub[i]);

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

      /************************************/
    case 2: /* Forward difference formula */
      /************************************/

      delta = MIN(deltay, deltap);
      rdelta = ONE/delta;
      r2delta = HALF/delta;
      
      /* Forward perturb y, y' and parameter */
      N_VLinearSum(delta, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(delta, ypsub[i], ONE, ypsub[0], yptemp);
      sdata->p[j-1] = psave + delta*pbar;

      /* Save residual in resvalsub[i] */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, resvalsub[i], rdata);
      if (retval != 0) return(retval);

      /* Estimate the residual for the i-th sensitivity equation */
      N_VLinearSum(rdelta, resvalsub[i], -rdelta, resvalsub[0], resvalsub[i]);

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

      /************************************************/
    case 3: /* Forward differences for separate terms */
      /************************************************/

      /* Forward perturb y and y' */
      N_VLinearSum(deltay, yysub[i], ONE, yysub[0], ytemp);
      N_VLinearSum(deltay, ypsub[i], ONE, ypsub[0], yptemp);

      /* Save residual in resvalsub[i] */
      retval = sdata->res_ptr(Ny, t, ytemp, yptemp, resvalsub[i], rdata);
      if (retval != 0) return(retval);

      /* Save the first difference quotient in resvalsub[i] */
      N_VLinearSum(rdeltay, resvalsub[i], -rdeltay, resvalsub[0], resvalsub[i]);

      /* Forward perturb parameter */
      sdata->p[j-1] = psave + deltap*pbar;

      /* Save residual in restemp */
      retval = sdata->res_ptr(Ny, t, yysub[0], ypsub[0], restemp, rdata);
      if (retval != 0) return(retval);

      /* Save the second difference quotient in restemp */
      N_VLinearSum(rdeltap, restemp, -rdeltap, resvalsub[0], restemp);

      /* Add the difference quotients for the residual for the i-th
         sensitivity */
      N_VLinearSum(ONE, resvalsub[i], ONE, restemp, resvalsub[i]);

      /* Restore original value of parameter */
      sdata->p[j-1] = psave;

      break;

    default:
      /* do nothing */
      break;
    }
  }
  return(SUCCESS);
}


/*******************************************************************/
/******** END Helper Function Implementation ***********************/
/*******************************************************************/


/***************************************************************/
/************** END SENSITIVITY Implementation *****************/
/***************************************************************/
