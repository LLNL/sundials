/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CPBANDPRE module, which
 * provides a banded difference quotient Jacobian-based
 * preconditioner and solver routines for use with CPSPGMR,
 * CPSPBCG, or CPSPTFQMR.
 *
 * Summary:
 * These routines provide a band matrix preconditioner based on
 * difference quotients of the ODE right-hand side function f
 * (for explicit-form ODEs) or on the residual function F (for
 * implicit-form ODEs).
 * The user supplies parameters
 *   mu = upper half-bandwidth (number of super-diagonals)
 *   ml = lower half-bandwidth (number of sub-diagonals)
 * The routines generate a band matrix of bandwidth ml + mu + 1
 * and use this to form a preconditioner for use with one of the
 * CSPILS iterative linear solvers. Although this matrix is intended
 * to approximate the Jacobian df/dy (respectively dF/dy + gamma*dF/dy,
 * it may be a very crude approximation. The true Jacobian need not 
 * be banded, or its true bandwith may be larger than ml + mu + 1, 
 * as long as the banded approximation generated here is sufficiently 
 * accurate to speed convergence as a preconditioner.
 *
 * Usage:
 *   The following is a summary of the usage of this module.
 *   Details of the calls to CPodeCreate, CPodeMalloc, CPSp*,
 *   and CPode are available in the User Guide.
 *   To use these routines, the sequence of calls in the user
 *   main program should be as follows:
 *
 *   #include <cpodes/cpodes_bandpre.h>
 *   #include <nvector_serial.h>
 *   ...
 *   void *bp_data;
 *   ...
 *   Set y0
 *   ...
 *   cpode_mem = CPodeCreate(...);
 *   ier = CPodeMalloc(...);
 *   ...
 *   flag = CPSptfqmr(cpode_mem, pretype, maxl);
 *     -or-
 *   flag = CPSpgmr(cpode_mem, pretype, maxl);
 *     -or-
 *   flag = CPSpbcg(cpode_mem, pretype, maxl);
 *   ...
 *   flag = CPBandPrecInit(cpode_mem, N, mu, ml);
 *   ...
 *   flag = CPode(...);
 *   ...
 *   Free y0
 *   ...
 *   CPodeFree(cpode_mem);
 *
 * Notes:
 * (1) Include this file for the CPBandPrecData type definition.
 * (2) In the CPBandPrecAlloc call, the arguments N is the 
 *     problem dimension.
 * (3) In the CPBPSp* call, the user is free to specify
 *     the input pretype and the optional input maxl.
 * -----------------------------------------------------------------
 */

#ifndef _CPBANDPRE_H
#define _CPBANDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CPBandPrecInit
 * -----------------------------------------------------------------
 * CPBandPrecInit allocates and initializes the BANDPRE preconditioner
 * module. This functino must be called AFTER one of the SPILS linear
 * solver modules has been attached to the CPODES integrator.
 *
 * The parameters of CPBandPrecInit are as follows:
 *
 * cpode_mem is the pointer to CPODES memory returned by CPodeCreate.
 *
 * N is the problem size.
 *
 * mu is the upper half bandwidth.
 *
 * ml is the lower half bandwidth.
 *
 * The return value of CPBandPrecInit is one of:
 *   CPSPILS_SUCCESS if no errors occurred
 *   CPSPILS_MEM_NULL if the integrator memory is NULL
 *   CPSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CPSPILS_ILL_INPUT if an input has an illegal value
 *   CPSPILS_MEM_FAIL if a memory allocation request failed
 *
 * NOTE: The band preconditioner assumes a serial implementation
 *       of the NVECTOR package. Therefore, CPBandPrecInit will
 *       first test for a compatible N_Vector internal
 *       representation by checking for required functions.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBandPrecInit(void *cpode_mem, int N, int mu, int ml);

/*
 * -----------------------------------------------------------------
 * Optional output functions : CPBandPrecGet*
 * -----------------------------------------------------------------
 * CPBandPrecGetWorkSpace returns the real and integer work space used
 *                        by CPBANDPRE.
 * CPBandPrecGetNumFctEvals returns the number of calls made from
 *                          CPBANDPRE to the user's ODE function.
 *
 * The return value of CPBandPrecGet* is one of:
 *    CPBANDPRE_SUCCESS    if successful
 *    CPBANDPRE_PDATA_NULL if the bp_data memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBandPrecGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CPBandPrecGetNumFctEvals(void *cpode_mem, long int *nfevalsBP);


#ifdef __cplusplus
}
#endif

#endif
