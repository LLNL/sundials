/******************************************************************
 * File          : fcvdense1.c                                    *
 * Programmers   : Radu Serban  and Alan Hindmarsh@ LLNL          *
 * Version of    : 30 March 2003                                  *
 *----------------------------------------------------------------*
 *                                                                *
 * Fortran/C interface routines for CVODE/CVDENSE, for the case   *
 * of a user-supplied Jacobian approximation routine.             *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvdense.h"       /* CVDense prototype, type DenseMat                */

/***************************************************************************/

/* Prototype of the Fortran routine */
void FCV_DJAC(integertype*, realtype*, realtype*, realtype*, realtype*, 
              realtype*, realtype*, realtype*);

/***************************************************************************/

void FCV_DENSE1(integertype *neq, int *ier)
{
  /* Call CVDense:
   *neq          is the problem dimension
     CV_cvodemem is the pointer to the CVODE memory block 
     CVDenseJac  is a pointer to the dense Jac routine
     NULL        is a pointer to jac_data                 */

  *ier = CVDense(CV_cvodemem, *neq, CVDenseJac, NULL);
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

void CVDenseJac(integertype N, DenseMat J, realtype t, 
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;

  ydata   = N_VGetData(y);
  fydata  = N_VGetData(fy);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  jacdata = DENSE_COL(J,0);

  FCV_DJAC(&N, &t, ydata, fydata, jacdata, v1data, v2data, v3data); 

}
