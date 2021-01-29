/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CPODES dense linear solver, CPDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CPDENSE_H
#define _CPDENSE_H

#include <cpodes/cpodes_direct.h>
#include <sundials/sundials_dense.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CPDense
 * -----------------------------------------------------------------
 * A call to the CPDense function links the main integrator with
 * the CPDENSE linear solver.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CPDense is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPDense(void *cpode_mem, int N);

/*
 * -----------------------------------------------------------------
 * Function : CPDenseProj
 * -----------------------------------------------------------------
 * A call to the CPDenseProj function links the main integrator with
 * the CPDENSE linear solver.
 *
 * cpode_mem  the pointer to the integrator memory returned by
 *            CPodeCreate.
 * Nc         the number of constraints
 * Ny         the number of states (size of the ODE system).
 * fact_type  the type of factorization used for the constraint
 *            Jcobian G. Legal values are CPDLS_LU and CPDLS_LQ.
 *
 * The return value of CPDense is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type);

#ifdef __cplusplus
}
#endif

#endif
