/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2008-04-18 19:42:35 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CPODES band linear solver, CPBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CPBAND_H
#define _CPBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_direct.h>
#include <sundials/sundials_band.h>

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
