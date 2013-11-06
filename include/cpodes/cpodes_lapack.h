/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2008-04-18 19:42:36 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * begincopyright(llns)
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * endcopyright(llns)
 * -----------------------------------------------------------------
 * Header file for the CPODES dense linear solver CPLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CPLAPACK_H
#define _CPLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_direct.h>
#include <sundials/sundials_lapack.h>

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CPLapackDense
 * -----------------------------------------------------------------
 * A call to the CPLapackDense function links the main integrator
 * with the CPLAPACK linear solver using dense Jacobians.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CPLapackDense is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPLapackDense(void *cpode_mem, int N);

/*
 * -----------------------------------------------------------------
 * Function : CPLapackBand
 * -----------------------------------------------------------------
 * A call to the CPLapackBand function links the main integrator
 * with the CPLAPACK linear solver using banded Jacobians. 
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian approximation.
 *
 * The return value of CPLapackBand is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing or
 *                       if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPLapackBand(void *cpode_mem, int N, int mupper, int mlower);

/*
 * -----------------------------------------------------------------
 * Function : CPLapackDenseProj
 * -----------------------------------------------------------------
 * A call to the CPLapackDenseProj function links the main integrator 
 * with the CPLAPACK linear solver using dense Jacobians.
 *
 * cpode_mem  the pointer to the integrator memory returned by
 *            CPodeCreate.
 * Nc         the number of constraints
 * Ny         the number of states (size of the ODE system).
 * fact_type  the type of factorization used for the constraint
 *            Jacobian G. Legal values are CPDLS_LU, CPDLS_QR,
 *            and CPDLS_SC.
 *
 * The return value of CPLapackDenseProj is one of:
 *    CPDLS_SUCCESS   if successful
 *    CPDLS_MEM_NULL  if the CPODES memory was NULL
 *    CPDLS_MEM_FAIL  if there was a memory allocation failure
 *    CPDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPLapackDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type);


#ifdef __cplusplus
}
#endif

#endif
