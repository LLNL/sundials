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
 * This is the header file for the CPODES scaled preconditioned GMRES 
 * linear solver, CPSPGMR.
 * -----------------------------------------------------------------
 */

#ifndef _CPSPGMR_H
#define _CPSPGMR_H

#include <cpodes/cpodes_spils.h>
#include <sundials/sundials_spgmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CPSpgmr
 * -----------------------------------------------------------------
 * A call to the CPSpgmr function links the main CPODES integrator
 * with the CPSPGMR linear solver.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * pretype   is the type of user preconditioning to be done.
 *           This must be one of the four enumeration constants
 *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined 
 *           in sundials_iterative.h.
 *           These correspond to no preconditioning,
 *           left preconditioning only, right preconditioning
 *           only, and both left and right preconditioning,
 *           respectively.
 *
 * maxl      is the maximum Krylov dimension. This is an
 *           optional input to the CPSPGMR solver. Pass 0 to
 *           use the default value CPSPGMR_MAXL=5.
 *
 * The return value of CPSpgmr is one of:
 *    CPSPILS_SUCCESS   if successful
 *    CPSPILS_MEM_NULL  if the CPODES memory was NULL
 *    CPSPILS_MEM_FAIL  if there was a memory allocation failure
 *    CPSPILS_ILL_INPUT if a required vector operation is missing
 * The above constants are defined in cpodes_spils.h
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPSpgmr(void *cpode_mem, int pretype, int maxl);


#ifdef __cplusplus
}
#endif

#endif
