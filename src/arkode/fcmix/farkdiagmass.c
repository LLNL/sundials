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
 * Fortran/C interface routines for ARKODE/ARKDENSE, for the case
 * of a user-supplied Jacobian approximation routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_direct.h>
#include <sunmatrix/sunmatrix_diagonal.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_DIAGMASS(realtype *T, realtype *DMASS, 
                            sunindextype *IPAR,  realtype *RPAR, 
                            realtype *V1, realtype *V2, realtype *V3, 
                            int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to ARKDlsSetMassFn; see 
   farkode.h for further details */
void FARK_DENSESETMASS(int *ier)
{
  *ier = ARKDlsSetMassFn(ARK_arkodemem, FARKDiagMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKDIAGMASS; see 
   farkode.h for additional information  */
int FARKDiagMass(realtype t, SUNMatrix M, void *user_data, 
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *massdata, *v1data, *v2data, *v3data;
  FARKUserData ARK_userdata;

  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  massdata = N_VGetArrayPointer(SUNDiagonalMatrix_Diag(M));
  ARK_userdata = (FARKUserData) user_data;

  FARK_DIAGMASS(&t, massdata, ARK_userdata->ipar, ARK_userdata->rpar, 
                v1data, v2data, v3data, &ier); 
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
