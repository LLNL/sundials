/******************************************************************
 *                                                                *
 * File          : sensfpvspgmr2.c                                *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * Fortran/C interface routines for sensitivity variant of        *
 * PVODE/CVSPGMR                                                  *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */
#include "sensfpvode.h" /* sensitivity version of Fortran/C interface      */
#include "senscvspgmr.h" /* SensCVSpgmr prototype                          */
#include "sensitivity.h" /* sensitivity data types and prototypes          */

/***************************************************************************/

void SFCV_SPGMR2 (int *pretype, int *gstype, int *maxl, real *delt)
{
  /* Call SensCVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     SCVPreco    is a pointer to the preconditioner setup routine
     SCVPSol     is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data                             */


  SensCVSpgmr (CV_cvodemem, *pretype, *gstype, *maxl, *delt, SCVPreco, 
	       SCVPSol, NULL);
}

/***************************************************************************/

/* C function SCVPreco to interface between CVODE and a Fortran subroutine
   PVPRECO for setup of a Krylov preconditioner.
   Addresses of Nlocal, t, jok, gamma, h, uround, y, fy, ewt, vtemp1, vtemp2, 
   vtemp3, the addresses jcurPtr and nfePtr, and p are passed to PVPRECO, 
   using the macros N_VLOCLENGTH and N_VDATA.  A return flag ier from PVPRECO
   is returned by SCVPreco.
   Auxiliary data is assumed to be communicated by Common. */

int SCVPreco(integer N, real t, N_Vector y, N_Vector fy, boole jok,
	     boole *jcurPtr, real gamma, N_Vector ewt, real h,
	     real uround, long int *nfePtr, void *P_data,
	     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  int ier = 0;
  integer Nlocal;
  CVodeMem cv_mem;
  SensData sdata;
  real *p;

  cv_mem = (CVodeMem) CV_cvodemem;
  sdata = (SensData) cv_mem->cv_f_data;
  p = sdata->p;

  Nlocal = N_VLOCLENGTH(y);
  ydata = N_VDATA(y);
  fydata = N_VDATA(fy);
  ewtdata = N_VDATA(ewt);
  v1data = N_VDATA(vtemp1);
  v2data = N_VDATA(vtemp2);
  v3data = N_VDATA(vtemp3);

  SFCV_PRECO (&Nlocal, &t, ydata, fydata, &jok, jcurPtr, &gamma, ewtdata,
           &h, &uround, nfePtr, v1data, v2data, v3data, &ier, p);

  return(ier);
}
