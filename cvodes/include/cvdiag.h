/*
 * -----------------------------------------------------------------
 * $Revision: 1.6.2.1 $
 * $Date: 2005-01-24 21:40:57 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES diagonal linear
 * solver, CVDIAG.
 * -----------------------------------------------------------------
 */

#ifndef _CVDIAG_H
#define _CVDIAG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Function : CVDiag
 * -----------------------------------------------------------------
 * A call to the CVDiag function links the main integrator with
 * the CVDIAG linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * The return value of CVDiag is one of:
 *    CVDIAG_SUCCESS   if successful
 *    CVDIAG_MEM_NULL  if the cvode memory was NULL
 *    CVDIAG_MEM_FAIL  if there was a memory allocation failure
 *    CVDIAG_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int CVDiag(void *cvode_mem);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVDIAG linear solver
 * -----------------------------------------------------------------
 *
 * CVDiagGetWorkSpace returns the real and integer workspace used
 *                    by CVDIAG.
 * CVDiagGetNumRhsEvals returns the number of calls to the user
 *                      f routine due to finite difference Jacobian
 *                      evaluation.
 *                      Note: The number of diagonal approximate
 *                      Jacobians formed is equal to the number of
 *                      CVDiagSetup calls. This number is available
 *                      through CVodeGetNumLinSolvSetups.
 * CVDiagGetLastFlag returns the last error flag set by any of
 *                   the CVDIAG interface functions.
 *
 * The return value of CVDiagGet* is one of:
 *    CVDIAG_SUCCESS   if successful
 *    CVDIAG_MEM_NULL  if the cvode memory was NULL
 *    CVDIAG_LMEM_NULL if the cvdiag memory was NULL
 * -----------------------------------------------------------------
 */

int CVDiagGetWorkSpace(void *cvode_mem, long int *lenrwDI, long int *leniwDI);
int CVDiagGetNumRhsEvals(void *cvode_mem, long int *nfevalsDI);
int CVDiagGetLastFlag(void *cvode_mem, int *flag);

#define CVDIAG_SUCCESS    0
#define CVDIAG_MEM_NULL  -1
#define CVDIAG_LMEM_NULL -2
#define CVDIAG_ILL_INPUT -3
#define CVDIAG_MEM_FAIL  -4

#define CVDIAG_INV_FAIL   1

#ifdef __cplusplus
}
#endif

#endif
