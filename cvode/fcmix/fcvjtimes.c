/*
 * -----------------------------------------------------------------
 * $Revision: 1.22 $
 * $Date: 2005-10-05 20:31:20 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * The C function FCVJtimes is to interface between the
 * CVSP* module and the user-supplied Jacobian-vector
 * product routine FCVJTIMES. Note the use of the generic name
 * FCV_JTIMES below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvsptfqmr.h"      /* CVSptfqmr prototype                            */
#include "cvspbcg.h"        /* CVSpbcg prototype                              */
#include "cvspgmr.h"        /* CVSpgmr prototype                              */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                               */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FCV_JTIMES(realtype*, realtype*, realtype*, realtype*, 
		       realtype*, realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_SPGMRSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSpgmrSetJacTimesVecFn(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVSpgmrSetJacTimesVecFn(CV_cvodemem, FCVJtimes, NULL);
  }
}

/***************************************************************************/

void FCV_SPBCGSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSpbcgSetJacTimesVecFn(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVSpbcgSetJacTimesVecFn(CV_cvodemem, FCVJtimes, NULL);
  }
}

/***************************************************************************/

void FCV_SPTFQMRSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSptfqmrSetJacTimesVecFn(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVSptfqmrSetJacTimesVecFn(CV_cvodemem, FCVJtimes, NULL);
  }
}

/***************************************************************************/

/* C function  FCVJtimes to interface between CVODE and  user-supplied
   Fortran routine FCVJTIMES for Jacobian * vector product.
   Addresses of v, Jv, t, y, fy, h, and work are passed to FCVJTIMES,
   using the routine N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVJTIMES is returned by FCVJtimes.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
              N_Vector y, N_Vector fy,
              void *jac_data, N_Vector work)
{
  realtype *vdata, *Jvdata, *ydata, *fydata, *wkdata;
  realtype h;

  int ier = 0;
  
  CVodeGetLastStep(CV_cvodemem, &h);

  vdata   = N_VGetArrayPointer(v);
  Jvdata  = N_VGetArrayPointer(Jv);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  wkdata  = N_VGetArrayPointer(work);

  FCV_JTIMES (vdata, Jvdata, &t, ydata, fydata, &h, wkdata, &ier);

  return(ier);
}
