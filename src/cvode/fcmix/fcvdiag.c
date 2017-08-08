/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVDLS, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_direct.h>
#include <sunmatrix/sunmatrix_diagonal.h>

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_DIAGJAC(realtype *T, realtype *Y, realtype *FY, 
                          realtype *DJAC, realtype *H, long int *IPAR, 
                          realtype *RPAR, realtype *V1, realtype *V2, 
                          realtype *V3, int *ier);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_DIAGSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVDlsSetJacFn(CV_cvodemem, NULL);
  } else {
    *ier = CVDlsSetJacFn(CV_cvodemem, FCVDiagJac);
  }
}

/***************************************************************************/

/* C function CVDiagJac interfaces between CVODE and a Fortran subroutine
   FCVDiagJAC for solution of a linear system with diagonal Jacobian 
   approximation.  Addresses of arguments are passed to FCVDIAGJAC, using 
   accessor functions from the SUNDiagonalMatrix and N_Vector modules. */

int FCVDiagJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector vtemp1, N_Vector vtemp2,
               N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  jacdata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(J));

  CV_userdata = (FCVUserData) user_data;

  FCV_DIAGJAC(&t, ydata, fydata, jacdata, &h, CV_userdata->ipar,
              CV_userdata->rpar, v1data, v2data, v3data, &ier); 
  return(ier);
}

