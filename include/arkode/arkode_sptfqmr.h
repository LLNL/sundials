/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the header file for the ARKODE scaled preconditioned 
 TFQMR linear solver, ARKSPTFQMR.
---------------------------------------------------------------*/

#ifndef _ARKSPTFQMR_H
#define _ARKSPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode.h>
#include <arkode/arkode_spils.h>
#include <sundials/sundials_sptfqmr.h>

/*---------------------------------------------------------------
 ARKSptfqmr:

 A call to the ARKSptfqmr function links the main ARKODE integrator
 with the ARKSPTFQMR linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 pretype   is the type of user preconditioning to be done.
           This must be one of the four enumeration constants
           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
           in iterative.h. These correspond to no preconditioning,
           left preconditioning only, right preconditioning
           only, and both left and right preconditioning,
           respectively.

 maxl      is the maximum Krylov dimension. This is an
           optional input to the ARKSPTFQMR solver. Pass 0 to
           use the default value ARKSPILS_MAXL=5.

 The return value of ARKSptfqmr is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_MEM_FAIL  if there was a memory allocation failure
    ARKSPILS_ILL_INPUT if a required vector operation is missing
 The above constants are defined in arkode_spils.h
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSptfqmr(void *arkode_mem, int pretype, int maxl);


/*---------------------------------------------------------------
 ARKMassSptfqmr:

 A call to the ARKMassSptfqmr function links the mass matrix solve
 with the ARKSPTFQMR linear solver module.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 pretype   is the type of user preconditioning to be done.
           This must be one of the four enumeration constants
           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
           in iterative.h. These correspond to no preconditioning,
           left preconditioning only, right preconditioning
           only, and both left and right preconditioning,
           respectively.

 maxl      is the maximum Krylov dimension. This is an
           optional input to the ARKSPTFQMR solver. Pass 0 to
           use the default value ARKSPILS_MAXL=5.

 mtimes    is the user-supplied mass-matrix-vector product 
           routine.

 The return value of ARKMassSptfqmr is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_MEM_FAIL  if there was a memory allocation failure
    ARKSPILS_ILL_INPUT if a required vector operation is missing
 The above constants are defined in arkode_spils.h
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassSptfqmr(void *arkode_mem, int pretype, 
				   int maxl, ARKMTimesFn mtimes,
				   void *mtimes_data);


#ifdef __cplusplus
}
#endif

#endif
