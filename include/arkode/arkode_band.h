/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Header file for the ARKODE band linear solver, ARKBAND.
---------------------------------------------------------------*/

#ifndef _ARKBAND_H
#define _ARKBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_direct.h>
#include <sundials/sundials_band.h>
 
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
