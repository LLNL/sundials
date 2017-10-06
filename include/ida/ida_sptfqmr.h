/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
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
 * This is the public header file for the IDA scaled preconditioned
 * TFQMR linear solver module, IDASPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPTFQMR_H
#define _IDASPTFQMR_H

#include <ida/ida_spils.h>
#include <sundials/sundials_sptfqmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmr
 * -----------------------------------------------------------------
 * A call to the IDASptfqmr function links the main integrator with
 * the IDASPTFQMR linear solver module. Its parameters are as
 * follows:
 *
 * IDA_mem  is the pointer to memory block returned by IDACreate.
 *
 * maxl     is the maximum Krylov subspace dimension, an
 *          optional input. Pass 0 to use the default value.
 *          Otherwise pass a positive integer.
 *
 * The return values of IDASptfqmr are:
 *    IDASPILS_SUCCESS    if successful
 *    IDASPILS_MEM_NULL   if the ida memory was NULL
 *    IDASPILS_MEM_FAIL   if there was a memory allocation failure
 *    IDASPILS_ILL_INPUT  if there was illegal input.
 * The above constants are defined in ida_spils.h
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASptfqmr(void *ida_mem, int maxl);


#ifdef __cplusplus
}
#endif

#endif
