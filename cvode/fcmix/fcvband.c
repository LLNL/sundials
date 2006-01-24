/*
 * -----------------------------------------------------------------
 * $Revision: 1.19 $
 * $Date: 2006-01-24 00:49:25 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVBAND, for the case of 
 * a user-supplied Jacobian approximation routine.                
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"            /* actual fn. names, prototypes and global vars.  */
#include "cvode_impl.h"        /* definition of CVodeMem type                    */
#include "cvode_band.h"        /* CVBand prototype                               */
#include "sundials_nvector.h"  /* definitions of type N_Vector and vector macros */
#include "sundials_types.h"    /* definition of type realtype                    */

/******************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_BJAC(long int*, long int*, long int*, long int*, /* N, MU, ML, EBAND */
                       realtype*, realtype*, realtype*,            /* T, Y, FY         */
                       realtype*,                                  /* BJAC             */
                       realtype*,                                  /* H                */
                       long int*, realtype*,                       /* IPAR, RPAR       */
                       realtype*, realtype*, realtype*,            /* V1, V2, V3       */
                       int*);                                      /* IER              */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_BANDSETJAC(int *flag, int *ier)
{
  CVodeMem cv_mem;

  if (*flag == 0) {
    *ier = CVBandSetJacFn(CV_cvodemem, NULL, NULL);
  } else {
    cv_mem = (CVodeMem) CV_cvodemem;
    *ier = CVBandSetJacFn(CV_cvodemem, FCVBandJac, cv_mem->cv_f_data);
  }
}

/***************************************************************************/

/* C function CVBandJac interfaces between CVODE and a Fortran subroutine
   FCVBJAC for solution of a linear system with band Jacobian approximation.
   Addresses of arguments are passed to FCVBJAC, using the macro 
   BAND_COL from BAND and the routine N_VGetArrayPointer from NVECTOR.
   The address passed for J is that of the element in column 0 with row 
   index -mupper.  An extended bandwith equal to (J->smu) + mlower + 1 is
   passed as the column dimension of the corresponding array.
   Auxiliary data is assumed to be communicated by Common. */

int FCVBandJac(long int N, long int mupper, long int mlower,
               BandMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int eband;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  eband = (J->smu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  CV_userdata = (FCVUserData) jac_data;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, jacdata, &h,
           CV_userdata->ipar, CV_userdata->rpar, v1data, v2data, v3data, &ier);

  return(0);
}
