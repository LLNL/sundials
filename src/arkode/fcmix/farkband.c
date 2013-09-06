/*---------------------------------------------------------------
  $Revision: 1.0 $
  $Date: $
 ---------------------------------------------------------------- 
  Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
  Fortran/C interface routines for ARKODE/ARKBAND, for the case 
  of a user-supplied Jacobian approximation routine.                
 --------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_band.h>


/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_BJAC(long int *N, long int *MU, long int *ML,
			long int *EBAND, realtype *T, 
			realtype *Y, realtype *FY, 
			realtype *BJAC, realtype *H,
			long int *IPAR, realtype *RPAR, 
			realtype *V1, realtype *V2, 
			realtype *V3, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to ARKDlsSetBandJacFn; see farkode.h 
   for further details */
void FARK_BANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKDlsSetBandJacFn(ARK_arkodemem, NULL);
  } else {
    *ier = ARKDlsSetBandJacFn(ARK_arkodemem, FARKBandJac);
  }
}

/*=============================================================*/

/* C interface to user-supplied Fortran subroutine FARKBJAC; see 
   farkode.h for further details */
int FARKBandJac(long int N, long int mupper, 
		long int mlower, realtype t, N_Vector y, 
		N_Vector fy, DlsMat J, void *user_data, 
		N_Vector vtemp1, N_Vector vtemp2, 
		N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int eband;
  FARKUserData ARK_userdata;

  ARKodeGetLastStep(ARK_arkodemem, &h);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  eband   = (J->s_mu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;
  ARK_userdata = (FARKUserData) user_data;

  FARK_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, 
	    jacdata, &h, ARK_userdata->ipar, ARK_userdata->rpar, 
	    v1data, v2data, v3data, &ier);

  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
