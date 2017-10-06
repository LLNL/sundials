/*
 * -----------------------------------------------------------------
 * $Revision: 4489 $
 * $Date: 2015-04-29 17:15:44 -0700 (Wed, 29 Apr 2015) $
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
 * Header file for the IDAS dense linear solver IDASLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _IDALAPACK_H
#define _IDALAPACK_H

#include <idas/idas_direct.h>
#include <sundials/sundials_lapack.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : IDALapackDense
 * -----------------------------------------------------------------
 * A call to the IDALapackDense function links the main integrator
 * with the IDALAPACK linear solver using dense Jacobians.
 *
 * ida_mem is the pointer to the integrator memory returned by
 *           IDACreate.
 *
 * N is the size of the DAE system.
 *
 * The return value of IDALapackDense is one of:
 *    IDADLS_SUCCESS   if successful
 *    IDADLS_MEM_NULL  if the IDAS memory was NULL
 *    IDADLS_MEM_FAIL  if there was a memory allocation failure
 *    IDADLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDALapackDense(void *ida_mem, int N);

/*
 * -----------------------------------------------------------------
 * Function : IDALapackBand
 * -----------------------------------------------------------------
 * A call to the IDALapackBand function links the main integrator
 * with the IDALAPACK linear solver using banded Jacobians. 
 *
 * ida_mem is the pointer to the integrator memory returned by
 *           IDACreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian approximation.
 *
 * The return value of IDALapackBand is one of:
 *    IDADLS_SUCCESS   if successful
 *    IDADLS_MEM_NULL  if the IDAS memory was NULL
 *    IDADLS_MEM_FAIL  if there was a memory allocation failure
 *    IDADLS_ILL_INPUT if a required vector operation is missing
 *                        or if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDALapackBand(void *ida_mem, int N, int mupper, int mlower);

/*
 * -----------------------------------------------------------------
 * Function: IDALapackDenseB
 * -----------------------------------------------------------------
 * IDALapackDenseB links the main IDAS integrator with the dense
 * IDALAPACK linear solver for the backward integration.
 * The 'which' argument is the int returned by IDACreateB.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDALapackDenseB(void *ida_mem, int which, int NeqB);

/*
 * -----------------------------------------------------------------
 * Function: IDALapackBandB
 * -----------------------------------------------------------------
 * IDALapackBandB links the main IDAS integrator with the band
 * IDALAPACK linear solver for the backward integration.
 * The 'which' argument is the int returned by IDACreateB.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDALapackBandB(void *ida_mem, int which, int NeqB, int mupperB, int mlowerB);


#ifdef __cplusplus
}
#endif

#endif
