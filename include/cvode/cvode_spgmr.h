/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-10-11 16:34:10 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE scaled preconditioned GMRES 
 * linear solver, CVSPGMR.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPGMR_H
#define _CVSPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvode/cvode_spils.h>
#include <sundials/sundials_spgmr.h>

/*
 * -----------------------------------------------------------------
 * Function : CVSpgmr
 * -----------------------------------------------------------------
 * A call to the CVSpgmr function links the main CVODE integrator
 * with the CVSPGMR linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
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
 *           optional input to the CVSPGMR solver. Pass 0 to
 *           use the default value CVSPGMR_MAXL=5.
 *
 * The return value of CVSpgmr is one of:
 *    CVSPILS_SUCCESS   if successful
 *    CVSPILS_MEM_NULL  if the cvode memory was NULL
 *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
 *    CVSPILS_ILL_INPUT if a required vector operation is missing
 * The above constants are defined in cvode_spils.h
 *
 * -----------------------------------------------------------------
 */

int CVSpgmr(void *cvode_mem, int pretype, int maxl);


#ifdef __cplusplus
}
#endif

#endif
