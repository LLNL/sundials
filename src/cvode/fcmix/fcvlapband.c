/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007-04-23 23:37:19 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVLAPACK, for the case
 * of a user-supplied band Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_lapack.h>

/***************************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_BJAC(int*, int*, int*, int*,           /* N, MU, ML, EBAND */
                       realtype*, realtype*, realtype*,  /* T, Y, FY         */
                       realtype*,                        /* LBJAC            */
                       realtype*,                        /* H                */
                       long int*, realtype*,             /* IPAR, RPAR       */
                       realtype*, realtype*, realtype*,  /* V1, V2, V3       */
                       int*);                            /* IER              */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_LAPACKBANDSETJAC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {

    *ier = CVDlsSetBandJacFn(CV_cvodemem, NULL);

  } else {

    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVDlsSetBandJacFn(CV_cvodemem, FCVLapackBandJac);

  }

}

/***************************************************************************/

/* The C function FCVLapackBandJac interfaces between CVODE and a 
 * Fortran subroutine FCVLBJAC for the solution of a linear system using
 * Lapack with band Jacobian approximation.
 * Addresses of arguments are passed to FCVLBJAC, using the macro 
 * LAPACK_BAND_COL and the routine N_VGetArrayPointer from NVECTOR.
 * The address passed for J is that of the element in column 0 with row 
 * index -mupper.  An extended bandwith equal to (J->smu) + mlower + 1 is
 * passed as the column dimension of the corresponding array.
 * Auxiliary data is assumed to be communicated by Common. 
 */

int FCVLapackBandJac(int N, int mupper, int mlower,
                     realtype t, N_Vector y, N_Vector fy, 
                     DlsMat J, void *f_data,
                     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  int eband;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  eband = (J->s_mu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  CV_userdata = (FCVUserData) f_data;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, jacdata, &h,
           CV_userdata->ipar, CV_userdata->rpar, v1data, v2data, v3data, &ier);

  return(ier);
}
