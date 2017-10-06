/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Header file for the KINSOL dense linear solver KINLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _KINLAPACK_H
#define _KINLAPACK_H

#include <kinsol/kinsol_direct.h>
#include <sundials/sundials_lapack.h>

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
 * Function : KINLapackDense
 * -----------------------------------------------------------------
 * A call to the KINLapackDense function links the main solver
 * with the KINLAPACK linear solver using dense Jacobians.
 *
 * kinmem is the pointer to the solver memory returned by KINCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of KINLapackDense is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINDLS_MEM_FAIL  if there was a memory allocation failure
 *    KINDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINLapackDense(void *kinmem, int N);

/*
 * -----------------------------------------------------------------
 * Function : KINLapackBand
 * -----------------------------------------------------------------
 * A call to the KINLapackBand function links the main solver
 * with the KINLAPACK linear solver using banded Jacobians. 
 *
 * kinmem is the pointer to the solver memory returned by KINCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian approximation.
 *
 * The return value of KINLapackBand is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINDLS_MEM_FAIL  if there was a memory allocation failure
 *    KINDLS_ILL_INPUT if a required vector operation is missing
 *                        or if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINLapackBand(void *kinmem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
