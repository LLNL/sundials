/*
 * -----------------------------------------------------------------
 * $Revision: 1.13 $
 * $Date: 2004-10-21 20:55:05 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * The C function FCVJtimes is to interface between the CVSPGMR module
 * and the user-supplied Jacobian-times-vector routine FCVJTIMES.     
 * Note the use of the generic name FCV_JTIMES below.                 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvspgmr.h"        /* CVSpgmr prototype                              */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                               */
#include "nvector.h"        /* definitions of type N_Vector and vector macros */
#include "sundialstypes.h"  /* definition of type realtype                    */

/* Prototype of the Fortran routine */
extern void FCV_JTIMES(realtype*, realtype*, realtype*, realtype*, 
		       realtype*, realtype*, realtype*, realtype*,
		       int*);

/***************************************************************************/

void FCV_SPGMRSETJAC(int *flag, int *ier)
{
  if (*flag == 0) CVSpgmrSetJacTimesVecFn(CV_cvodemem, NULL);
  else CVSpgmrSetJacTimesVecFn(CV_cvodemem, FCVJtimes);
}

/***************************************************************************/

/* C function  FCVJtimes to interface between CVODE and  user-supplied
   Fortran routine FCVJTIMES for Jacobian * vector product.
   Addresses of v, Jv, t, y, fy, vnrm, ewt, h, uround, ytemp, and
   the address nfePtr, are passed to FCVJTIMES, using the routine
   N_VGetArrayPointer from NVECTOR. A return flag ier from FCVJTIMES is 
   returned by FCVJtimes.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
              N_Vector y, N_Vector fy,
              void *jac_data, N_Vector work)
{
  N_Vector ewt;
  realtype *vdata, *Jvdata, *ydata, *fydata, *ewtdata, *wkdata;
  realtype h;

  int ier = 0;

  CVodeGetErrWeights(CV_cvodemem, &ewt);
  CVodeGetLastStep(CV_cvodemem, &h);

  vdata   = N_VGetArrayPointer(v);
  Jvdata  = N_VGetArrayPointer(Jv);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  wkdata  = N_VGetArrayPointer(work);
  ewtdata = N_VGetArrayPointer(ewt);

  FCV_JTIMES (vdata, Jvdata, &t, ydata, fydata, 
              ewtdata, &h, wkdata, &ier);

  return(ier);
}
