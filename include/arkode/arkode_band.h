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
 * Header file for the ARKODE band linear solver, ARKBAND.
 *--------------------------------------------------------------*/

#ifndef _ARKBAND_H
#define _ARKBAND_H

#include <arkode/arkode_direct.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
 ARKBand:

 A call to the ARKBand function links the main ARKODE integrator
 with the ARKBAND linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 mupper is the upper bandwidth of the band Jacobian
        approximation.

 mlower is the lower bandwidth of the band Jacobian
        approximation.

 The return value of ARKBand is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the arkode memory was NULL
    ARKDLS_MEM_FAIL  if there was a memory allocation failure
    ARKDLS_ILL_INPUT if a required vector operation is missing or
                       if a bandwidth has an illegal value.
---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKBand(void *arkode_mem, long int N, 
			    long int mupper, long int mlower);

/*---------------------------------------------------------------
 ARKMassBand:

 A call to the ARKMassBand function links the mass matrix solve
 with the ARKBAND linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 mupper is the upper bandwidth of the band Jacobian
        approximation.

 mlower is the lower bandwidth of the band Jacobian
        approximation.

 bmass is the user-supplied band mass matrix setup routine.

 The return value of ARKBand is one of:
    ARKDLS_SUCCESS   if successful
    ARKDLS_MEM_NULL  if the arkode memory was NULL
    ARKDLS_MEM_FAIL  if there was a memory allocation failure
    ARKDLS_ILL_INPUT if a required vector operation is missing or
                       if a bandwidth has an illegal value.
---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKMassBand(void *arkode_mem, long int N, 
				long int mupper, long int mlower,
				ARKDlsBandMassFn bmass);

#ifdef __cplusplus
}
#endif

#endif
