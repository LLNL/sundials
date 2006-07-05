/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:27:50 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE scaled preconditioned 
 * Bi-CGSTAB linear solver, CVSPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPBCG_H
#define _CVSPBCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvode/cvode_spils.h>
#include <sundials/sundials_spbcgs.h>

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcg
 * -----------------------------------------------------------------
 * A call to the CVSpbcg function links the main CVODE integrator
 * with the CVSPBCG linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
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
 *           optional input to the CVSPBCG solver. Pass 0 to
 *           use the default value CVSPBCG_MAXL=5.
 *
 * The return value of CVSpbcg is one of:
 *    CVSPILS_SUCCESS   if successful
 *    CVSPILS_MEM_NULL  if the cvode memory was NULL
 *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
 *    CVSPILS_ILL_INPUT if a required vector operation is missing
 * The above constants are defined in cvode_spils.h
 *
 * -----------------------------------------------------------------
 */

int CVSpbcg(void *cvode_mem, int pretype, int maxl);


#ifdef __cplusplus
}
#endif

#endif
