/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * The C functions FARKMTSetup and FARKMtimes are to interface 
 * between the ARKLS and ARKLSMASS modules and the user-supplied 
 * mass-matrix-vector setup/product routines FARKMTSETUP and 
 * FARKJTIMES. Note the use of the generic names FARK_MTSETUP 
 * and FARK_MTIMES in the code below.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_MTSETUP(realtype *T, long int *IPAR, 
                           realtype *RPAR, int *IER);
  extern void FARK_MTIMES(realtype *V, realtype *MV, realtype *T, 
                          long int *IPAR, realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetMassTimes; see 
   farkode.h for further information */
void FARK_SPILSSETMASS(int *ier)
{ FARK_LSSETMASS(ier); }

/* Fortran interface to C routine ARKStepSetMassTimes; see 
   farkode.h for further information */
void FARK_LSSETMASS(int *ier)
{
  ARKodeMem ark_mem;
  ark_mem = (ARKodeMem) ARK_arkodemem;
  *ier = ARKStepSetMassTimes(ARK_arkodemem, FARKMTSetup, 
                             FARKMtimes, ark_mem->user_data);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKMTSETUP; see
   farkode.h for further information */
int FARKMTSetup(realtype t, void *user_data)
{
  FARKUserData ARK_userdata;
  int ier = 0;
  ARK_userdata = (FARKUserData) user_data;
  FARK_MTSETUP(&t, ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}

/* C interface to user-supplied Fortran routine FARKMTIMES; see
   farkode.h for further information */
int FARKMtimes(N_Vector v, N_Vector Mv, realtype t, void *user_data)
{
  realtype *vdata, *Mvdata;
  FARKUserData ARK_userdata;
  int ier = 0;
  
  vdata  = N_VGetArrayPointer(v);
  Mvdata = N_VGetArrayPointer(Mv);
  ARK_userdata = (FARKUserData) user_data;
  FARK_MTIMES(vdata, Mvdata, &t, ARK_userdata->ipar, 
              ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
