/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2008-04-18 19:42:36 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE band linear solver, CVBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CVBAND_H
#define _CVBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvode/cvode_direct.h>
#include <sundials/sundials_band.h>
 
/*
 * -----------------------------------------------------------------
 * Function : CVBand
 * -----------------------------------------------------------------
 * A call to the CVBand function links the main CVODE integrator
 * with the CVBAND linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian
 *        approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian
 *        approximation.
 *
 * The return value of CVBand is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the cvode memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing or
 *                       if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBand(void *cvode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
