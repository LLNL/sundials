/******************************************************************
 * File          : fcvband1.c                                     *
 * Programmers   : Radu Serban and Alan Hindmarsh @ LLNL          *
 * Version of    : 22 July 2002                                   *
 *----------------------------------------------------------------*
 *                                                                *
 * Fortran/C interface routines for CVODE/CVBAND, for the case of *
 * a user-supplied Jacobian approximation routine.                *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype  */
#include "nvector.h"       /* definitions of type N_Vector and vector macros */
#include "fcvode.h"        /* actual function names, prototypes, global vars.*/
#include "cvband.h"        /* CVBand prototype                               */

/******************************************************************************/

/* Prototype of the Fortran routine */
void FCV_BJAC(integertype*, integertype*, integertype*, integertype*, realtype*,
              realtype*, realtype*, realtype*, realtype*, realtype*, realtype*,
              long int*,  realtype*, realtype*, realtype*);

/******************************************************************************/

void FCV_BAND1(integertype *mupper, integertype *mlower, int *ier)
{
  /* Call CVBand:
     CV_cvodemem is the pointer to the CVODE memory block 
     *mupper     is the upper bandwidth
     *mlower     is the lower bandwidth
     CVBandJac   is a pointer to the band Jac routine
     NULL        is a pointer to jac_data                 */

  *ier = CVBand(CV_cvodemem, *mupper, *mlower, CVBandJac, NULL);
}

/***************************************************************************/

void FCV_REINBAND1(integertype *mupper, integertype *mlower, int *ier)
{
  /* Call CVReInitBand:
     CV_cvodemem is the pointer to the CVODE memory block 
     *mupper     is the upper bandwidth
     *mlower     is the lower bandwidth
     CVBandJac   is a pointer to the band Jac routine
     NULL        is a pointer to jac_data                 */

  *ier = CVReInitBand(CV_cvodemem, *mupper, *mlower, CVBandJac, NULL);
}

/***************************************************************************/

/* C function CVBandJac interfaces between CVODE and a Fortran subroutine
   CVBJAC for solution of a linear system with band Jacobian approximation.
   Addresses of arguments are passed to CVBJAC, using the macro 
   BAND_COL from BAND and the routine N_VGetData from NVECTOR.
   The address passed for J is that of the element in column 0 with row 
   index -mupper.  An extended bandwith equal to (J->smu) + mlower + 1 is
   passed as the column dimension of the corresponding array.
   Auxiliary data is assumed to be communicated by Common. */

void CVBandJac(integertype N, integertype mupper, integertype mlower,
               BandMat J, RhsFn f, void *f_data, realtype t,
               N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
               realtype uround, void *jac_data, long int *nfePtr,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *ewtdata, *jacdata, *v1data, *v2data, *v3data;
  integertype eband;

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  eband = (J->smu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, ewtdata, 
           &h, &uround, jacdata, nfePtr, v1data, v2data, v3data);

}
