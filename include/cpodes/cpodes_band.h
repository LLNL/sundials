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
 * This is the header file for the CPODES band linear solver, CPBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CPBAND_H
#define _CPBAND_H

#include <cpodes/cpodes_direct.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CPBand
 * -----------------------------------------------------------------
 * A call to the CPBand function links the main CPODES integrator
 * with the CPBAND linear solver.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian
 *        approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian
 *        approximation.
 *
 * The return value of CPBand is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing
 *                       or if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPBand(void *cpode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
