/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CPODES scaled preconditioned TFQMR 
 * linear solver, CPSPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _CPSPTFQMR_H
#define _CPSPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_spils.h>
#include <sundials/sundials_sptfqmr.h>

/*
 * -----------------------------------------------------------------
 * Function : CPSptfqmr
 * -----------------------------------------------------------------
 * A call to the CPSptfqmr function links the main CPODES integrator
 * with the CPSPTFQMR linear solver.
 *
 * cpode_mem is the pointer to the integrator memory returned by
 *           CPodeCreate.
 *
 * pretype   is the type of user preconditioning to be done.
 *           This must be one of the four enumeration constants
 *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
 *           in iterative.h. These correspond to no preconditioning,
 *           left preconditioning only, right preconditioning
 *           only, and both left and right preconditioning,
 *           respectively.
 *
 * maxl      is the maximum Krylov dimension. This is an
 *           optional input to the CPSPTFQMR solver. Pass 0 to
 *           use the default value CPSPILS_MAXL=5.
 *
 * The return value of CPSptfqmr is one of:
 *    CPSPILS_SUCCESS   if successful
 *    CPSPILS_MEM_NULL  if the CPODES memory was NULL
 *    CPSPILS_MEM_FAIL  if there was a memory allocation failure
 *    CPSPILS_ILL_INPUT if a required vector operation is missing
 * The above constants are defined in cpodes_spils.h
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CPSptfqmr(void *cpode_mem, int pretype, int maxl);

#ifdef __cplusplus
}
#endif

#endif
