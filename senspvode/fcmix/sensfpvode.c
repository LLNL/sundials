/******************************************************************
 *                                                                *
 * File          : sensfpvode.c                                   *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * Fortran/C interface routines for SensPVODE                     *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "cvode.h"    /* CVODE constants and prototypes                    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "mpi.h"      /* MPI types and constants                           */
#include "sensfpvode.h"  /* sensitivity version of Fortran/C interface     */
#include "sensitivity.h" /* sensitivity data types and prototypes          */

/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* SFPVMalloc/SFPVReInit Error Messages */

#define SFPVM         "SFPVMalloc/SFPVReInit-- "

#define MSG_REI_NO_MEM   "SensCVReInit-- cvode_mem = NULL illegal.\n\n"
#define MSG_REI_NO_SDATA "SensCVReInit-- s_data = NULL illegal.\n\n"
#define MSG_BAD_N        SFPVM "N=%ld < 1 illegal.\n\n"
#define MSG_BAD_NY       SFPVM "Ny=%ld < 1 illegal.\n\n"
#define MSG_BAD_NS       SFPVM "Ns=%ld < 0 illegal.\n\n"
#define MSG_P_NULL       SFPVM "p=NULL illegal.\n\n"
#define MSG_PBAR_NULL    SFPVM "pbar=NULL illegal.\n\n"
#define MSG_F_NULL       SFPVM "f=NULL illegal.\n\n"
#define MSG_FTEMP_NULL   SFPVM "ftemp=NULL illegal.\n\n"
#define MSG_YTEMP_NULL   SFPVM "ytemp=NULL illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/

/***************************************************************************/

void SFPV_MALLOC (integer *ny, integer *ns, integer *ntotal, real *t0, 
		  real *y0, integer *meth, integer *itmeth, integer *iatol, 
		  real *rtol, real *atol, integer *optin, long int *iopt, 
		  real *ropt, int *ier, real *p, real *pbar, real *rhomax)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;
  pData psave;

  CV_yvec = SensN_VMAKE(*ntotal, y0, PV_machEnv);
  if (CV_yvec == NULL) {
    *ier = -1;
    return;
  }
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { 
      atolvec = SensN_VMAKE(*ntotal, atol, PV_machEnv);
      if (atolvec == NULL) {
	*ier = -1;
	return;
      }
      itol = SV; atolptr = atolvec; 
    }

  psave = (pData) malloc(sizeof(*psave));
  psave->p = p;

  /* Call SensCVodeMalloc to initialize CVODE: 
     *ny     is the number of state variables
     *ns     is the number of sensitivity vectors
     *ntotal is the total problem size
     SensCVf is the right-hand side function in y'=f(t,y,p)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector of length *ntotal
     lmm     is the method specifier
     iter    is the iteration method specifier
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     psave   is the pointer to a data block that points to p (used in SensCVf)
     NULL    is the pointer to the error message file
     *optin  is the flag indicating optional inputs in iopt and ropt
     iopt    is an array used to communicate optional integer input and output
     ropt    is an array used to communicate optional real input and output
     PV_machEnv is the pointer to the machine environment block
     p       is a pointer to the array of real parameters in y' = f(t,y,p)
     pbar    is a pointer to an array used to scale the sensitivity vectors
     rhomax  is a real value used for selecting a finite difference formula
               for estimating the scaled sensitivity derivatives

     A pointer to CVODE problem memory is returned and stored in CV_cvodemem. */

  CV_cvodemem = SensCVodeMalloc(*ny, *ns, *ntotal, SensCVf, *t0, CV_yvec,
                         lmm, iter, itol, rtol, atolptr, psave, NULL,
                         *optin, iopt, ropt, PV_machEnv, p, pbar, *rhomax);

  *ier = (CV_cvodemem == NULL) ? -1 : 0 ;

}

/***************************************************************************/

void SFPV_REINIT (real *t0, real *y0, integer *meth, integer *itmeth,
     integer *iatol, real *rtol, real *atol, integer *optin,
     long int *iopt, real *ropt, int *ier, real *p, real *pbar, real *rhomax)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;  
  CVodeMem cv_mem;
  SensData sdata; 
  pData psave;
  integer Ntotal;

  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }
  cv_mem = (CVodeMem) CV_cvodemem;

  if (CV_yvec == NULL) {
    *ier = -1;
    return;
  }
  else {
    N_VDATA(CV_yvec) = y0;
  }
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { 
      Ntotal = cv_mem->cv_N;
      atolvec = SensN_VMAKE(Ntotal, atol, PV_machEnv);
      if (atolvec == NULL) {
	*ier = -1;
	return;
      }
      itol = SV; atolptr = atolvec; 
    }

  /* Call SensCVReInit to re-initialize CVODE: 
     CV_cvodemem is the pointer to the CVODE memory block 
     SensCVf is the user's right-hand side function in y'=f(t,y,p)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector of length *ntotal
     lmm     is the method specifier
     iter    is the iteration method specifier
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     psave   is the pointer to the data block that points to p (used in SensCVf)
     NULL    is the pointer to the error message file
     *optin  is the flag indicating optional inputs in iopt and ropt
     iopt    is an array used to communicate optional integer input and output
     ropt    is an array used to communicate optional real input and output
     PV_machEnv is the pointer to the machine environment block
     p       is a pointer to the array of real parameters in y' = f(t,y,p)
     pbar    is a pointer to an array used to scale the sensitivity vectors
     rhomax  is a real value used for selecting a finite difference formula
               for estimating the scaled sensitivity derivatives

     A pointer to CVODE problem memory is returned and stored in CV_cvodemem. */

  *ier = SensCVReInit(CV_cvodemem, SensCVf, *t0, CV_yvec, lmm, iter, 
                 itol, rtol, atolptr, psave, NULL, *optin, iopt, ropt, 
                 PV_machEnv, p, pbar, *rhomax);

  /* psave and its parameter array will be accessed by SensCVf.           */    

  sdata = (SensData) cv_mem->cv_f_data;
  psave = (pData) sdata->f_data;
  psave->p = p;
  
}

/***************************************************************************/

void SFCV_SPGMR0 (int *gstype, int *maxl, real *delt)
{
  /* Call SensCVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     0           is the preconditioner type (none)
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     NULL        is a pointer to the preconditioner setup routine
     NULL        is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data                             */

  SensCVSpgmr (CV_cvodemem, 0, *gstype, *maxl, *delt, NULL, NULL, NULL);
}

/***************************************************************************/

void SFCV_FREE ()
{
  CVodeMem cv_mem;
  SensData sdata;

  cv_mem = (CVodeMem) CV_cvodemem;
  sdata = (SensData) cv_mem->cv_f_data;
  
  free(sdata->f_data);

  /* Call SensCVodeFree:
     CV_cvodemem is the pointer to the CVODE memory block */

  SensCVodeFree (CV_cvodemem);
}

/***************************************************************************/

/* C function SensCVf to interface between CVODE and a Fortran subroutine 
   PVFUN.
   Addresses of Nlocal, t, y, ydot, and p are passed to PVFUN, 
   using the  macros N_VLOCLENGTH and N_VDATA from the NVECTOR module.
   Auxiliary data is assumed to be communicated by Common. */

void SensCVf(integer N, real t, N_Vector y, N_Vector ydot, void *p_save)
{
  real *ydata, *dydata;
  integer Nlocal;
  real *p;
  pData psave;

  psave = (pData) p_save;
  p = psave->p;

  Nlocal = N_VLOCLENGTH(y);
  ydata = N_VDATA(y);
  dydata = N_VDATA(ydot);

  SFCV_FUN (&Nlocal, &t, ydata, dydata, p);
}
