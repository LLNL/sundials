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
 * The C function FARKPSet is to interface between the ARKSP*
 * modules and the user-supplied preconditioner setup routine 
 * FARKPSET. Note the use of the generic name FARK_PSET in the 
 * code below.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_spils.h>

/*=============================================================*/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_MASSPSET(realtype *T, long int *IPAR, 
			    realtype *RPAR, realtype *W1, 
			    realtype *W2, realtype *W3, int *IER);
  extern void FARK_MASSPSOL(realtype *T, realtype *R, realtype *Z, 
			    realtype *DELTA, int *LR, long int *IPAR, 
			    realtype *RPAR, realtype *WRK, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKSpilsSetMassPreconditioner; see 
   farkode.h for further details */
void FARK_SPILSSETMASSPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKSpilsSetMassPreconditioner(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKSpilsSetMassPreconditioner(ARK_arkodemem, 
					 FARKMassPSet, FARKMassPSol);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKMASSPSET; see 
   farkode.h for further details */
int FARKMassPSet(realtype t, void *user_data, N_Vector vtemp1, 
		 N_Vector vtemp2, N_Vector vtemp3)
{
  int ier = 0;
  realtype *v1data, *v2data, *v3data;
  FARKUserData ARK_userdata;

  v1data = N_VGetArrayPointer(vtemp1);
  v2data = N_VGetArrayPointer(vtemp2);
  v3data = N_VGetArrayPointer(vtemp3);
  ARK_userdata = (FARKUserData) user_data;

  FARK_MASSPSET(&t, ARK_userdata->ipar, ARK_userdata->rpar,
		v1data, v2data, v3data, &ier);
  return(ier);
}


/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKMASSPSOL; see 
   farkode.h for further details */
int FARKMassPSol(realtype t, N_Vector r, N_Vector z, realtype delta,
		 int lr, void *user_data, N_Vector vtemp)
{
  int ier = 0;
  realtype *vtdata, *rdata, *zdata;
  FARKUserData ARK_userdata;

  vtdata = N_VGetArrayPointer(vtemp);
  rdata  = N_VGetArrayPointer(r);
  zdata  = N_VGetArrayPointer(z);
  ARK_userdata = (FARKUserData) user_data;

  FARK_MASSPSOL(&t, rdata, zdata, &delta, &lr, ARK_userdata->ipar, 
		ARK_userdata->rpar, vtdata, &ier);
  return(ier);
}


/*===============================================================
   EOF
===============================================================*/
