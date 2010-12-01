/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2010-12-01 22:27:37 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVDENSE, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_dense.h>

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_DJAC(long int*,                        /* N          */
                       realtype*, realtype*, realtype*,  /* T, Y, FY   */
                       realtype*,                        /* DJAC       */
                       realtype*,                        /* H          */ 
                       long int*, realtype*,             /* IPAR, RPAR */
                       realtype*, realtype*, realtype*,  /* V1, V2, V3 */
                       int *ier);                        /* IER        */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_DENSESETJAC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {
    *ier = CVDlsSetDenseJacFn(CV_cvodemem, NULL);
  } else {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVDlsSetDenseJacFn(CV_cvodemem, FCVDenseJac);
  }
}

/***************************************************************************/

/* C function CVDenseJac interfaces between CVODE and a Fortran subroutine
   FCVDJAC for solution of a linear system with dense Jacobian approximation.
   Addresses of arguments are passed to FCVDJAC, using the macro 
   DENSE_COL from DENSE and the routine N_VGetArrayPointer from NVECTOR.
   Auxiliary data is assumed to be communicated by Common. */

int FCVDenseJac(long int N, realtype t, 
                N_Vector y, N_Vector fy, 
                DlsMat J, void *user_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
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

  jacdata = DENSE_COL(J,0);

  CV_userdata = (FCVUserData) user_data;

  FCV_DJAC(&N, &t, ydata, fydata, jacdata, &h, 
           CV_userdata->ipar, CV_userdata->rpar, v1data, v2data, v3data, &ier); 

  return(ier);
}

