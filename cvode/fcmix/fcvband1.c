/******************************************************************
 * File          : fcvband1.c                                     *
 * Programmers   : Radu Serban @ LLNL                             *
 * Version of    : 28 March 2002                                  *
 *----------------------------------------------------------------*
 *                                                                *
 * Fortran/C interface routine for CVODE/CVBAND, for the case of  *
 * user supplied Jacobian approximation routine.                  *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fcvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvband.h"   /* CVBand prototype                                 */


/***************************************************************************/

/* Prototypes of the Fortran routines */
void FCV_BJAC(integer*, integer*, integer*, real*, real*, real*, real*, 
              real*, real*, real*, long int*, real*, real*, real*);

/***************************************************************************/

void FCV_BAND1(integer *mupper, integer *mlower, int *ier)
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

void FCV_REINBAND1(integer *mupper, integer *mlower, int *ier)
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
   Auxiliary data is assumed to be communicated by Common. */

void CVBandJac(integer N, integer mupper, integer mlower,
               BandMat J, RhsFn f, void *f_data, real t,
               N_Vector y, N_Vector fy, N_Vector ewt, real h,
               real uround, void *jac_data, long int *nfePtr,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  real *jacdata;

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  jacdata = BAND_COL(J,0);

  FCV_BJAC(&N, &mupper, &mlower, jacdata, &t, ydata, fydata, 
           ewtdata, &h, &uround, nfePtr, v1data, v2data, v3data);

}
