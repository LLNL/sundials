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
 * Fortran/C interface routines for ARKODE/ARKLAPACKBAND, for the 
 * case of a user-supplied Jacobian approximation routine.                
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_lapack.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_BMASS(long int *N, long int *MU, long int *ML,
			 long int *EBAND, realtype *T, 
			 realtype *BMASS, long int *IPAR, 
			 realtype *RPAR, realtype *V1, 
			 realtype *V2, realtype *V3, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to ARKDlsSetBandMassFn; see farkode.h 
   for further details */
void FARK_LAPACKBANDSETMASS(int *ier)
{
  *ier = ARKDlsSetBandMassFn(ARK_arkodemem, FARKLapackBandMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran subroutine FARKBMASS; see 
   farkode.h for further details */
int FARKLapackBandMass(long int N, long int mupper, 
		       long int mlower, realtype t, DlsMat M, 
		       void *user_data, N_Vector vtemp1, 
		       N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *massdata, *v1data, *v2data, *v3data;
  long int eband;
  FARKUserData ARK_userdata;

  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  eband   = (M->s_mu) + mlower + 1;
  massdata = BAND_COL(M,0) - mupper;
  ARK_userdata = (FARKUserData) user_data;

  FARK_BMASS(&N, &mupper, &mlower, &eband, &t, massdata, 
	     ARK_userdata->ipar, ARK_userdata->rpar, v1data, 
	     v2data, v3data, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
