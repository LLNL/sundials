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
 * Fortran/C interface routines for CVODE/CVDENSE, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvdense.h"        /* CVDense prototype and type DenseMat            */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                               */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */



/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_DJAC(long int*, realtype*, realtype*, realtype*, realtype*, 
                       realtype*, realtype*, realtype*, realtype*);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVDenseSetJacFn(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVDenseSetJacFn(CV_cvodemem, FCVDenseJac, NULL);
  }
}

/***************************************************************************/

/* C function CVDenseJac interfaces between CVODE and a Fortran subroutine
   FCVDJAC for solution of a linear system with dense Jacobian approximation.
   Addresses of arguments are passed to FCVDJAC, using the macro 
   DENSE_COL from DENSE and the routine N_VGetArrayPointer from NVECTOR.
   Auxiliary data is assumed to be communicated by Common. */

void FCVDenseJac(long int N, DenseMat J, realtype t, 
                 N_Vector y, N_Vector fy, void *jac_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  jacdata = DENSE_COL(J,0);

  FCV_DJAC(&N, &t, ydata, fydata, jacdata, &h, v1data, v2data, v3data); 

}

