/******************************************************************
 * File          : fcvdense1.c                                    *
 * Programmers   : Radu Serban @ LLNL                             *
 * Version of    : 28 March 2002                                  *
 *----------------------------------------------------------------*
 *                                                                *
 * Fortran/C interface routine for CVODE/CVDENSE, for the case of *
 * user supplied Jacobian approximation routine.                  *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fcvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvdense.h"  /* CVDense prototype                                 */

/***************************************************************************/

/* Prototypes of the Fortran routines */
void FCV_DJAC(integer*, real*, real*, real*, real*, real*, real*, real*,
              long int*, real*, real*, real*);

/***************************************************************************/

void FCV_DENSE1(int *ier)
{
  /* Call CVDense:
     CV_cvodemem is the pointer to the CVODE memory block 
     CVDenseJac  is a pointer to the dense Jac routine
     NULL        is a pointer to jac_data                 */

  *ier = CVDense(CV_cvodemem, CVDenseJac, NULL);
}

/***************************************************************************/

void FCV_REINDENSE1(int *ier)
{
  /* Call CVReInitDense:
     CV_cvodemem is the pointer to the CVODE memory block 
     CVDenseJac  is a pointer to the dense Jac routine
     NULL        is a pointer to jac_data                 */

  *ier = CVReInitDense(CV_cvodemem, CVDenseJac, NULL);
}

/***************************************************************************/

/* C function CVDenseJac interfaces between CVODE and a Fortran subroutine
   CVDJAC for solution of a linear system with dense Jacobian approximation.
   Addresses of arguments are passed to CVDJAC, using the macro 
   DENSE_COL from DENSE and the routine N_VGetData from NVECTOR.
   Auxiliary data is assumed to be communicated by Common. */

void CVDenseJac(integer N, DenseMat J, RhsFn f, void *f_data,
                real t, N_Vector y, N_Vector fy, N_Vector ewt,
                real h, real uround, void *jac_data,
                long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  real *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  real *jacdata;

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  jacdata = DENSE_COL(J,0);

  FCV_DJAC(&N, jacdata, &t, ydata, fydata, ewtdata, &h, &uround, 
           nfePtr, v1data, v2data, v3data); 

}
