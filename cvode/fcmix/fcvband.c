/******************************************************************
 * File          : fcvband.c                                     *
 * Programmers   : Radu Serban and Alan Hindmarsh @ LLNL          *
 * Version of    : 1 August 2003                                  *
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
void FCV_BJAC(integertype*, integertype*, integertype*, integertype*, 
              realtype*, realtype*, realtype*, realtype*,
              realtype*, realtype*, realtype*);


/***************************************************************************/

void FCV_BANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) CVBandSetJacFn(CV_cvodemem, NULL);
  else            CVBandSetJacFn(CV_cvodemem, CVBandJac);
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
               BandMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  integertype eband;

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  eband = (J->smu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  FCV_BJAC(&N, &mupper, &mlower, &eband, 
           &t, ydata, fydata, jacdata, 
           v1data, v2data, v3data);

}
