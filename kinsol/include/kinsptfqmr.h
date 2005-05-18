/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-05-18 18:17:42 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPTFQMR linear solver module header file (public version)
 * -----------------------------------------------------------------
 */

#ifndef _KINSPTFQMR_H
#define _KINSPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "sptfqmr.h"
#include "kinspils.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmr
 * -----------------------------------------------------------------
 * KINSptfqmr links the main KINSOL solver module with the SPTFQMR
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSptfqmrInit, KINSptfqmrSetup, KINSptfqmrSolve,
 * and KINSptfqmrFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPTFQMR_MAXL (predefined constant) to be used)
 *
 * If successful, KINSptfqmr returns KINSPTFQMR_SUCCESS. If an error
 * occurs, then KINSptfqmr returns an error code (negative integer
 * value).
 * -----------------------------------------------------------------
 */

int KINSptfqmr(void *kinmem, int maxl);

/*
 * -----------------------------------------------------------------
 * KINSptfqmr Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSptfqmr subroutine are the
 * following:
 *
 * KINSPTFQMR_SUCCESS : means the KINSPTFQMR linear solver module
 *                      (implementation of the TFQMR method) was
 *                      successfully initialized - allocated system
 *                      memory and set shared variables to default
 *                      values [0]
 *
 * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer
 *                       was given (must call the KINCreate and
 *                       KINMalloc memory allocation subroutines
 *                       prior to calling KINSptfqmr) [-1]
 *
 * KINSPTFQMR_MEM_FAIL : means either insufficient system resources
 *                       were available to allocate memory for the
 *                       main KINSPTFQMR data structure (type
 *                       KINSptfqmrMemRec), or the SptfqmrMalloc
 *                       subroutine failed (unable to allocate enough
 *                       system memory for vector storate and/or the
 *                       main SPTFQMR data structure
 *                       (type SptfqmrMemRec)) [-4]
 *
 * KINSPTFQMR_ILL_INPUT : means either a supplied parameter was invalid,
 *                        or the NVECTOR implementation is NOT
 *                        compatible [-3]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions (KINSPTFQMR)
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *       Function Name         |  Optional Input  [Default Value]
 *                             |
 * -----------------------------------------------------------------
 *                             |
 * KINSptfqmrSetPreconditioner | used to set the following:
 *                             |   (a) name of user-supplied routine
 *                             |       used to compute a preconditioner
 *                             |       matrix for the given linear
 *                             |       system (pset)
 *                             |       [NULL]
 *                             |   (b) name of user-supplied routine
 *                             |       used to apply preconditioner to
 *                             |       linear system (psolve)
 *                             |       [NULL]
 *                             |   (c) pointer to user-allocated system
 *                             |       memory that is passed to the pset
 *                             |       and psolve routines
 *                             |       [NULL]
 *                             |
 * KINSptfqmrSetJacTimesVecFn  | used to set the following:
 *                             |   (a) name of user-supplied subroutine
 *                             |       used to compute the matrix-vector
 *                             |       product J(u)*v, where J denotes
 *                             |       the system Jacobian (jtimes)
 *                             |       [KINSptfqmrDQJtimes] (see kinsptfqmr.c)
 *                             |   (b) pointer to a user-allocated memory
 *                             |       block that is passed to the jtimes
 *                             |       routine
 *                             |       [NULL]
 * -----------------------------------------------------------------
 */

int KINSptfqmrSetPreconditioner(void *kinmem,
				KINSpilsPrecSetupFn pset,
				KINSpilsPrecSolveFn psolve,
                                void *P_data);
int KINSptfqmrSetJacTimesVecFn(void *kinmem,
			       KINSpilsJacTimesVecFn jtimes,
			       void *J_data);

/*
 * -----------------------------------------------------------------
 * KINSptfqmrSet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSptfqmrSet* subroutines
 * are the following:
 *
 * KINSPTFQMR_SUCCESS : means the associated variable was successfully
 *                      set [0]
 *
 * KINSPTFQMR_ILL_INPUT : means the supplied parameter was invalid
 *                        (check error message) [-3]
 *
 * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                       given [-1]
 *
 * KINSPTFQMR_LMEM_NULL : means system memory has not yet been allocated
 *                        for SPTFQMR (lmem == NULL) [-2]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions (KINSPTFQMR)
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSPTFQMR linear
 * solver:
 *
 *        Function Name        |      Returned Value
 *                             |
 * -----------------------------------------------------------------
 *                             |
 * KINSptfqmrGetWorkSpace      | returns both integer workspace size
 *                             | (total number of long int-sized blocks
 *                             | of memory allocated by KINSPTFQMR for
 *                             | vector storage), and real workspace
 *                             | size (total number of realtype-sized
 *                             | blocks of memory allocated by KINSPTFQMR
 *                             | for vector storage)
 *                             |
 * KINSptfqmrGetNumPrecEvals   | total number of preconditioner
 *                             | evaluations (number of calls made
 *                             | to the user-defined pset routine)
 *                             |
 * KINSptfqmrGetNumPrecSolves  | total number of times preconditioner
 *                             | was applied to linear system (number
 *                             | of calls made to the user-supplied
 *                             | psolve function)
 *                             |
 * KINSptfqmrGetNumLinIters    | total number of linear iterations
 *                             | performed
 *                             |
 * KINSptfqmrGetNumConvFails   | total number of linear convergence
 *                             | failures
 *                             |
 * KINSptfqmrGetNumJtimesEvals | total number of times the matrix-
 *                             | vector product J(u)*v was computed
 *                             | (number of calls made to the jtimes
 *                             | subroutine)
 *                             |
 * KINSptfqmrGetNumFuncEvals   | total number of evaluations of the
 *                             | system function F(u) (number of
 *                             | calls made to the user-supplied
 *                             | func routine by the KINSPTFQMR module
 *                             | member subroutines)
 *                             |
 * KINSptfqmrGetLastFlag       | returns last flag returned by the
 *                             | linear solver
 * -----------------------------------------------------------------
 */

int KINSptfqmrGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG);
int KINSptfqmrGetNumPrecEvals(void *kinmem, long int *npevals);
int KINSptfqmrGetNumPrecSolves(void *kinmem, long int *npsolves);
int KINSptfqmrGetNumLinIters(void *kinmem, long int *nliters);
int KINSptfqmrGetNumConvFails(void *kinmem, long int *nlcfails);
int KINSptfqmrGetNumJtimesEvals(void *kinmem, long int *njvevals);
int KINSptfqmrGetNumFuncEvals(void *kinmem, long int *nfevalsSG); 
int KINSptfqmrGetLastFlag(void *kinmem, int *flag);

/*
 * -----------------------------------------------------------------
 * KINSptfqmrGet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSptfqmrGet* subroutines
 * are the following:
 *
 * KINSPTFQMR_SUCCESS : means the routine exited normally [0]
 *
 * KINSPTFQMR_ILL_INPUT : means at least one input parameter was
 *                        invalid (check error message(s)) [-3]
 *
 * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                       given [-1]
 *
 * KINSPTFQMR_LMEM : means a NULL KINSPTFQMR memory block pointer was
 *                   given [-2]
 * -----------------------------------------------------------------
 */

#define KINSPTFQMR_SUCCESS 0

#define KINSPTFQMR_MEM_NULL  -1
#define KINSPTFQMR_LMEM_NULL -2
#define KINSPTFQMR_ILL_INPUT -3
#define KINSPTFQMR_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
