/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * The C function FARKJtimes is to interface between the ARKSP* 
 * modules and the user-supplied Jacobian-vector product routine
 * FARKJTIMES. Note the use of the generic name FARK_JTIMES in
 * the code below.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_spils.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_JTIMES(realtype *V, realtype *JV, realtype *T, 
			  realtype *Y, realtype *FY, realtype *H,
			  long int *IPAR, realtype *RPAR,
			  realtype *WRK, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKSpilsSetJacTimesVecFn; see 
   farkode.h for further information */
void FARK_SPILSSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKSpilsSetJacTimesVecFn(ARK_arkodemem, NULL);
  } else {
    *ier = ARKSpilsSetJacTimesVecFn(ARK_arkodemem, FARKJtimes);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKJTIMES; see
   farkode.h for further information */
int FARKJtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
	       N_Vector fy, void *user_data, N_Vector work)
{
  realtype *vdata, *Jvdata, *ydata, *fydata, *wkdata;
  realtype h;
  FARKUserData ARK_userdata;
  int ier = 0;
  
  ARKodeGetLastStep(ARK_arkodemem, &h);

  vdata  = N_VGetArrayPointer(v);
  Jvdata = N_VGetArrayPointer(Jv);
  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  wkdata = N_VGetArrayPointer(work);

  ARK_userdata = (FARKUserData) user_data;
 
  FARK_JTIMES(vdata, Jvdata, &t, ydata, fydata, &h, ARK_userdata->ipar, 
	      ARK_userdata->rpar, wkdata, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
