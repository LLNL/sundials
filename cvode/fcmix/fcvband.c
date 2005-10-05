/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2005-10-05 20:31:20 $
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

#include "cvband.h"         /* CVBand prototype                               */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                               */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/******************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_BJAC(long int*, long int*, long int*, long int*, 
                       realtype*, realtype*, realtype*, realtype*,
                       realtype*, realtype*, realtype*, realtype*);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_BANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVBandSetJacFn(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVBandSetJacFn(CV_cvodemem, FCVBandJac, NULL);
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

void FCVBandJac(long int N, long int mupper, long int mlower,
                BandMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int eband;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  eband = (J->smu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, jacdata, 
           &h, v1data, v2data, v3data);

}
