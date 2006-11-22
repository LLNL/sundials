/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-11-22 00:12:47 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the IDAS dense linear solver IDASLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _IDALAPACK_H
#define _IDALAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <idas/idas_direct.h>
#include <sundials/sundials_lapack.h>

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
   * N is the size of the ODE system.
   *
   * The return value of IDALapackDense is one of:
   *    IDADIRECT_SUCCESS   if successful
   *    IDADIRECT_MEM_NULL  if the IDAS memory was NULL
   *    IDADIRECT_MEM_FAIL  if there was a memory allocation failure
   *    IDADIRECT_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int IDALapackDense(void *ida_mem, int N);

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
   *    IDADIRECT_SUCCESS   if successful
   *    IDADIRECT_MEM_NULL  if the IDAS memory was NULL
   *    IDADIRECT_MEM_FAIL  if there was a memory allocation failure
   *    IDADIRECT_ILL_INPUT if a required vector operation is missing
   *                        or if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int IDALapackBand(void *ida_mem, int N, int mupper, int mlower);

  /*
   * -----------------------------------------------------------------
   * Function: IDALapackDenseB
   * -----------------------------------------------------------------
   * IDALapackDenseB links the main IDAS integrator with the dense
   * IDALAPACK linear solver for the backward integration.
   * -----------------------------------------------------------------
   */

  int IDALapackDenseB(void *idaadj_mem, int NeqB);

  /*
   * -----------------------------------------------------------------
   * Function: IDALapackBandB
   * -----------------------------------------------------------------
   * IDALapackBandB links the main IDAS integrator with the band
   * IDALAPACK linear solver for the backward integration.
   * -----------------------------------------------------------------
   */

  int IDALapackBandB(void *idaadj_mem, int NeqB, int mupperB, int mlowerB);


#ifdef __cplusplus
}
#endif

#endif
