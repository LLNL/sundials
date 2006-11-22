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
 * Header file for the CVODE dense linear solver CVLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CVLAPACK_H
#define _CVLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvode/cvode_direct.h>
#include <sundials/sundials_lapack.h>

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CVLapackDense
   * -----------------------------------------------------------------
   * A call to the CVLapackDense function links the main integrator
   * with the CVLAPACK linear solver using dense Jacobians.
   *
   * cvode_mem is the pointer to the integrator memory returned by
   *           CVodeCreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of CVLapackDense is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
   *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CVLAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CVLapackDense(void *cvode_mem, int N);

  /*
   * -----------------------------------------------------------------
   * Function : CVLapackBand
   * -----------------------------------------------------------------
   * A call to the CVLapackBand function links the main integrator
   * with the CVLAPACK linear solver using banded Jacobians. 
   *
   * cvode_mem is the pointer to the integrator memory returned by
   *           CVodeCreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian approximation.
   *
   * The return value of CVLapackBand is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
   *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CVLAPACK_ILL_INPUT if a required vector operation is missing or
   *                       if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
