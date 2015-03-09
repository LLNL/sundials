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
 * This is the header file for the ARKODE scaled preconditioned 
 * GMRES linear solver, ARKSPGMR.
 *--------------------------------------------------------------*/

#ifndef _ARKSPGMR_H
#define _ARKSPGMR_H

#include <arkode/arkode.h>
#include <arkode/arkode_spils.h>
#include <sundials/sundials_spgmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 ARKSpgmr:

 A call to the ARKSpgmr function links the main ARKODE integrator
 with the ARKSPGMR linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 pretype   is the type of user preconditioning to be done.
           This must be one of the four enumeration constants
           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined 
           in sundials_iterative.h.
           These correspond to no preconditioning,
           left preconditioning only, right preconditioning
           only, and both left and right preconditioning,
           respectively.

 maxl      is the maximum Krylov dimension. This is an
           optional input to the ARKSPGMR solver. Pass 0 to
           use the default value ARKSPGMR_MAXL=5.

 The return value of ARKSpgmr is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_MEM_FAIL  if there was a memory allocation failure
    ARKSPILS_ILL_INPUT if a required vector operation is missing
 The above constants are defined in arkode_spils.h

---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSpgmr(void *arkode_mem, int pretype, int maxl);


/*---------------------------------------------------------------
 ARKMassSpgmr:

 A call to the ARKMassSpgmr function links the mass matrix solve
 with the ARKSPGMR linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 pretype   is the type of user preconditioning to be done.
           This must be one of the four enumeration constants
           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined 
           in sundials_iterative.h.
           These correspond to no preconditioning,
           left preconditioning only, right preconditioning
           only, and both left and right preconditioning,
           respectively.

 maxl      is the maximum Krylov dimension. This is an
           optional input to the ARKSPGMR solver. Pass 0 to
           use the default value ARKSPGMR_MAXL=5.

 mtimes    is the user-supplied mass-matrix-vector product 
           routine.

 The return value of ARKMassSpgmr is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_MEM_FAIL  if there was a memory allocation failure
    ARKSPILS_ILL_INPUT if a required vector operation is missing
 The above constants are defined in arkode_spils.h

---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassSpgmr(void *arkode_mem, int pretype, int maxl, 
				 ARKSpilsMassTimesVecFn mtimes,
				 void *mtimes_data);


#ifdef __cplusplus
}
#endif

#endif
