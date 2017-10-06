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
 * Header file for the ARKODE dense linear solver ARKLAPACK.
 *--------------------------------------------------------------*/

#ifndef _ARKLAPACK_H
#define _ARKLAPACK_H

#include <arkode/arkode_direct.h>
#include <sundials/sundials_lapack.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  EXPORTED FUNCTIONS
===============================================================*/

/*---------------------------------------------------------------
 ARKLapackDense:

 A call to the ARKLapackDense function links the main integrator
 with the ARKLAPACK linear solver using dense Jacobians.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 The return value of ARKLapackDense is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKLapackDense(void *arkode_mem, int N);


/*---------------------------------------------------------------
 ARKMassLapackDense:

 A call to the ARKMassLapackDense function links the mass matrix 
 solve with the ARKLAPACK linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 dmass is the user-supplied dense mass matrix setup routine.

 The return value of ARKMassLapackDense is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassLapackDense(void *arkode_mem, int N, 
				       ARKDlsDenseMassFn dmass);


/*---------------------------------------------------------------
 ARKLapackBand:

 A call to the ARKLapackBand function links the main integrator
 with the ARKLAPACK linear solver using banded Jacobians. 

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 mupper is the upper bandwidth of the band Jacobian approximation.

 mlower is the lower bandwidth of the band Jacobian approximation.

 The return value of ARKLapackBand is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing 
                        or if a bandwidth has an illegal value.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKLapackBand(void *arkode_mem, int N, int mupper, int mlower);

/*---------------------------------------------------------------
 ARKMassLapackBand:

 A call to the ARKMassLapackBand function links the mass matrix 
 solve with the ARKLAPACK linear solver. 

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 mupper is the upper bandwidth of the band Jacobian approximation.

 mlower is the lower bandwidth of the band Jacobian approximation.

 bmass is the user-supplied band mass matrix setup routine.

 The return value of ARKLapackBand is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing 
                        or if a bandwidth has an illegal value.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassLapackBand(void *arkode_mem, int N, int mupper, 
				      int mlower, ARKDlsBandMassFn bmass);

#ifdef __cplusplus
}
#endif

#endif
