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
 * Header file for the ARKODE dense linear solver, ARKDENSE.
 *---------------------------------------------------------------*/

#ifndef _ARKDENSE_H
#define _ARKDENSE_H

#include <arkode/arkode_direct.h>
#include <sundials/sundials_dense.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 Function: ARKDense
-----------------------------------------------------------------
 A call to the ARKDense function links the main integrator with
 the ARKDENSE linear solver.

 arkode_mem is the pointer to the integrator memory returned by
            ARKodeCreate.

 N is the size of the ODE system.

 The return value of ARKDense is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the arkode memory was NULL
    ARKDLS_MEM_FAIL  if there was a memory allocation failure
    ARKDLS_ILL_INPUT if a required vector operation is missing
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKDense(void *arkode_mem, long int N);

/*---------------------------------------------------------------
 Function: ARKMassDense
-----------------------------------------------------------------
 A call to the ARKMassDense function links the mass matrix solve
 with the ARKDENSE linear solver.

 arkode_mem is the pointer to the integrator memory returned by
            ARKodeCreate.

 N is the size of the ODE system.

 dmass is the user-supplied dense mass matrix setup function.

 The return value of ARKMassDense is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the arkode memory was NULL
    ARKDLS_MEM_FAIL  if there was a memory allocation failure
    ARKDLS_ILL_INPUT if a required vector operation is missing
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassDense(void *arkode_mem, long int N, 
				 ARKDlsDenseMassFn dmass);

#ifdef __cplusplus
}
#endif

#endif
