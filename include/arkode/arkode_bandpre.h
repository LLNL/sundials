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
 * This is the header file for the ARKBANDPRE module, which
 * provides a banded difference quotient Jacobian-based
 * preconditioner and solver routines for use with ARKSPGMR,
 * ARKSPBCG, ARKSPTFQMR, ARKSPFGMR or ARKPCG.
 *
 * Summary:
 * These routines provide a band matrix preconditioner based on
 * difference quotients of the implicit portion of the ODE 
 * right-hand side function f.  The user supplies parameters
 *   mu = upper half-bandwidth (number of super-diagonals)
 *   ml = lower half-bandwidth (number of sub-diagonals)
 * The routines generate a band matrix of bandwidth ml + mu + 1
 * and use this to form a preconditioner for use with the Krylov
 * linear solver in ARKSP*. Although this matrix is intended to
 * approximate the Jacobian dfi/dy, it may be a very crude
 * approximation. The true Jacobian need not be banded, or its
 * true bandwith may be larger than ml + mu + 1, as long as the
 * banded approximation generated here is sufficiently accurate
 * to speed convergence as a preconditioner.
 *
 * Usage:
 *   The following is a summary of the usage of this module.
 *   Details of the calls to ARKodeCreate, ARKSp*,
 *   and ARKode are available in the User Guide.
 *   To use these routines, the sequence of calls in the user
 *   main program should be as follows:
 *
 *   #include <arkode/arkode_bandpre.h>
 *   #include <nvector_serial.h>
 *   ...
 *   Set y0
 *   ...
 *   arkode_mem = ARKodeCreate();
 *   ier = ARKodeInit(...);
 *   ...
 *   flag = ARKSptfqmr(arkode_mem, pretype, maxl);
 *     -or-
 *   flag = ARKSpgmr(arkode_mem, pretype, maxl);
 *     -or-
 *   flag = ARKSpbcg(arkode_mem, pretype, maxl);
 *     -or-
 *   flag = ARKSpfgmr(arkode_mem, pretype, maxl);
 *     -or-
 *   flag = ARKPcg(arkode_mem, pretype, maxl);
 *   ...
 *   flag = ARKBandPrecInit(arkode_mem, N, mu, ml);
 *   ...
 *   flag = ARKode(...);
 *   ...
 *   Free y0
 *   ...
 *   ARKodeFree(&arkode_mem);
 *
 * Notes:
 * (1) Include this file for the ARKBandPrecData type definition.
 * (2) In the ARKBandPrecAlloc call, the arguments N is the
 *     problem dimension.
 * (3) In the ARKBPSp* call, the user is free to specify
 *     the input pretype and the optional input maxl.
 *--------------------------------------------------------------*/

#ifndef _ARKBANDPRE_H
#define _ARKBANDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 ARKBandPrecInit:

 ARKBandPrecInit allocates and initializes the BANDPRE preconditioner
 module. This functino must be called AFTER one of the SPILS linear
 solver modules has been attached to the ARKODE integrator.

 The parameters of ARKBandPrecInit are as follows:

 arkode_mem is the pointer to ARKODE memory returned by ARKodeCreate.

 N is the problem size.

 mu is the upper half bandwidth.

 ml is the lower half bandwidth.

 The return value of ARKBandPrecInit is one of:
   ARKSPILS_SUCCESS if no errors occurred
   ARKSPILS_MEM_NULL if the integrator memory is NULL
   ARKSPILS_LMEM_NULL if the linear solver memory is NULL
   ARKSPILS_ILL_INPUT if an input has an illegal value
   ARKSPILS_MEM_FAIL if a memory allocation request failed

 NOTE: The band preconditioner assumes a serial implementation
       of the NVECTOR package. Therefore, ARKBandPrecInit will
       first test for a compatible N_Vector internal
       representation by checking for required functions.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKBandPrecInit(void *arkode_mem, long int N, 
				    long int mu, long int ml);

/*---------------------------------------------------------------
 Optional output functions : ARKBandPrecGet*

 ARKBandPrecGetWorkSpace returns the real and integer work space used
                        by ARKBANDPRE.

 ARKBandPrecGetNumRhsEvals returns the number of calls made from
                          ARKBANDPRE to the user's right-hand side
                          routine f.

 The return value of ARKBandPrecGet* is one of:
   ARKSPILS_SUCCESS if no errors occurred
   ARKSPILS_MEM_NULL if the integrator memory is NULL
   ARKSPILS_LMEM_NULL if the linear solver memory is NULL
   ARKSPILS_PMEM_NULL if the preconditioner memory is NULL
---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKBandPrecGetWorkSpace(void *arkode_mem, 
					    long int *lenrwLS, 
					    long int *leniwLS);
SUNDIALS_EXPORT int ARKBandPrecGetNumRhsEvals(void *arkode_mem, 
					      long int *nfevalsBP);


#ifdef __cplusplus
}
#endif

#endif
