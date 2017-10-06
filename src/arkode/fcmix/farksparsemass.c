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
 * Fortran/C interface routines for ARKODE/ARKKLU, for the case
 * of a user-supplied sparse Jacobian routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_sparse.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_SPMASS(realtype *T, int *N, int *NNZ, 
			  realtype *MDATA, int *MRVALS, 
			  int *MCPTRS, long int *IPAR, 
			  realtype *RPAR, realtype *V1, 
			  realtype *V2, realtype *V3, int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKSlsSetSparseMassFn; see 
   farkode.h for further information */
void FARK_SPARSESETMASS(int *ier)
{
  *ier = ARKSlsSetSparseMassFn(ARK_arkodemem, FARKSparseMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKSPMASS; see 
   farkode.h for additional information  */
int FARKSparseMass(realtype t, SlsMat MassMat, void *user_data, 
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *v1data, *v2data, *v3data;
  FARKUserData ARK_userdata;

  v1data = N_VGetArrayPointer(vtemp1);
  v2data = N_VGetArrayPointer(vtemp2);
  v3data = N_VGetArrayPointer(vtemp3);
  ARK_userdata = (FARKUserData) user_data;

  FARK_SPMASS(&t, &(MassMat->NP), &(MassMat->NNZ), MassMat->data, 
	      MassMat->indexvals, MassMat->indexptrs, ARK_userdata->ipar, 
	      ARK_userdata->rpar, v1data, v2data, v3data, &ier); 
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
