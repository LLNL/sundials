/******************************************************************
 *                                                                *
 * File          : sensfpvspgmr1.c                                *
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

void SFCV_SPGMR1 (int *pretype, int *gstype, int *maxl, real *delt)
{
  /* Call SensCVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     NULL        is a pointer to the preconditioner setup routine
     SCVPSol     is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data                             */

  SensCVSpgmr (CV_cvodemem, *pretype, *gstype, *maxl, *delt, NULL, 
	       SCVPSol, NULL);
}

/***************************************************************************/

/* C function SCVPSol to interface between CVODE and a Fortran subroutine
   PVPSOL for solution of a Krylov preconditioner.
   Addresses of Nlocal, t, gamma, delta, lr, y, fy, vtemp, ewt, r, z, the
   address nfePtr, and p are passed to PVPSOL, using the macros N_VLOCLENGTH
   and N_VDATA.  A return flag ier from PVPSOL is returned by SCVPSol.
   Auxiliary data is assumed to be communicated by Common. */

int SCVPSol(integer N, real t, N_Vector y, N_Vector fy, N_Vector vtemp,
           real gamma, N_Vector ewt, real delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z)
{
  real *ydata, *fydata, *vtdata, *ewtdata, *rdata, *zdata;
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
  vtdata = N_VDATA(vtemp);
  ewtdata = N_VDATA(ewt);
  rdata = N_VDATA(r);
  zdata = N_VDATA(z);

  SFCV_PSOL (&Nlocal, &t, ydata, fydata, vtdata, &gamma, ewtdata, &delta,
	     nfePtr, rdata, &lr, zdata, &ier, p);

  return(ier);
}
