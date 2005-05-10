/*
 * -----------------------------------------------------------------
 * $Revision: 1.22 $
 * $Date: 2005-05-10 23:19:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSPGMR linear solver module header file
 * -----------------------------------------------------------------
 */

#ifndef _KINSPGMR_H
#define _KINSPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "spgmr.h"
#include "kinspils.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmr
 * -----------------------------------------------------------------
 * KINSpgmr links the main KINSOL solver module with the SPGMR
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSpgmrInit, KINSpgmrSetup, KINSpgmrSolve,
 * and KINSpgmrFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPGMR_MAXL (predefined constant) to be used)
 *
 * KINSpgmr return values: KINSPGMR_SUCCESS, KINSPGMR_MEM_NULL,
 * KINSPGMR_MEM_FAIL and KINSPGMR_ILL_INPUT (see below).
 * -----------------------------------------------------------------
 */

int KINSpgmr(void *kinmem, int maxl);

/*
 * -----------------------------------------------------------------
 * KINSpgmr Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpgmr subroutine are the
 * following:
 *
 * KINSPGMR_SUCCESS : means the KINSPGMR linear solver module
 *                    (implementation of the GMRES method) was
 *                    successfully initialized - allocated system
 *                    memory and set shared variables to default
 *                    values [0]
 *
 * KINSPGMR_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                     given (must call the KINCreate and KINMalloc
 *                     memory allocation subroutines prior to
 *                     calling KINSpgmr) [-1]
 *
 * KINSPGMR_MEM_FAIL : means either insufficient system resources
 *                     were available to allocate memory for the main
 *                     KINSPGMR data structure (type KINSpgmrMemRec),
 *                     or the SpgmrMalloc subroutine failed (unable
 *                     to allocate enough system memory for vector
 *                     storage and/or the main SPGMR data structure
 *                     (type SpgmrMemRec)) [-4]
 *
 * KINSPGMR_ILL_INPUT : means a supplied parameter was invalid
 *                      (check error message) [-3]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions (KINSPGMR)
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *       Function Name       |   Optional Input  [Default Value]
 *                           |
 * -----------------------------------------------------------------
 *                           |
 * KINSpgmrSetMaxRestarts    | maximum number of times the SPGMR
 *                           | (scaled preconditioned GMRES) linear
 *                           | solver can be restarted
 *                           | [0]
 *                           |
 * KINSpgmrSetPreconditioner | used to set the following:
 *                           |   (a) name of user-supplied routine
 *                           |       used to compute a preconditioner
 *                           |       matrix for the given linear
 *                           |       system (pset)
 *                           |       [NULL]
 *                           |   (b) name of user-supplied routine
 *                           |       used to apply preconditioner to
 *                           |       linear system (psolve)
 *                           |       [NULL]
 *                           |   (c) pointer to user-allocated system
 *                           |       memory that is passed to the pset
 *                           |       and psolve routines
 *                           |       [NULL]
 *                           |
 * KINSpgmrSetJacTimesVecFn  | used to set the following:
 *                           |   (a) name of user-supplied subroutine
 *                           |       used to compute the matrix-vector
 *                           |       product J(u)*v, where J denotes
 *                           |       the system Jacobian (jtimes)
 *                           |       [KINSpgmrDQJtimes] (see kinspgmr.c)
 *                           |   (b) pointer to a user-allocated memory
 *                           |       block that is passed to the jtimes
 *                           |       routine
 *                           |       [NULL]
 * -----------------------------------------------------------------
 */

int KINSpgmrSetMaxRestarts(void *kinmem, int maxrs);
int KINSpgmrSetPreconditioner(void *kinmem,
			      KINSpilsPrecSetupFn pset,
			      KINSpilsPrecSolveFn psolve,
			      void *P_data);
int KINSpgmrSetJacTimesVecFn(void *kinmem,
			     KINSpilsJacTimesVecFn jtimes,
			     void *J_data);

/*
 * -----------------------------------------------------------------
 * KINSpgmrSet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpgmrSet* subroutines
 * are the following:
 *
 * KINSPGMR_SUCCESS : means the associated parameter was successfully
 *                    set [0]
 *
 * KINSPGMR_ILL_INPUT : means the supplied parameter was invalid
 *                      (check error message) [-3]
 *
 * KINSPGMR_MEM_NULL : means a NULL KINSOL memory block pointer
 *                     was given [-1]
 *
 * KINSPGMR_LMEM_NULL : means system memory has not yet been
 *                      allocated for SPGMR (lmem == NULL) [-2]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions (KINSPGMR)
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSPGMR linear
 * solver:
 *
 *        Function Name       |      Returned Value
 *                            |
 * -----------------------------------------------------------------
 *                            |
 * KINSpgmrGetWorkSpace       | returns both integer workspace size
 *                            | (total number of long int-sized blocks
 *                            | of memory allocated by KINSPGMR for
 *                            | vector storage), and real workspace
 *                            | size (total number of realtype-sized
 *                            | blocks of memory allocated by KINSPGMR
 *                            | for vector storage)
 *                            |
 * KINSpgmrGetNumPrecEvals    | total number of preconditioner
 *                            | evaluations (number of calls made
 *                            | to the user-defined pset routine)
 *                            |
 * KINSpgmrGetNumPrecSolves   | total number of times preconditioner
 *                            | was applied to linear system (number
 *                            | of calls made to the user-supplied
 *                            | psolve function)
 *                            |
 * KINSpgmrGetNumLinIters     | total number of linear iterations
 *                            | performed
 *                            |
 * KINSpgmrGetNumConvFails    | total number of linear convergence
 *                            | failures
 *                            |
 * KINSpgmrGetNumJtimesEvals  | total number of times the matrix-
 *                            | vector product J(u)*v was computed
 *                            | (number of calls made to the jtimes
 *                            | subroutine)
 *                            |
 * KINSpgmrGetNumFuncEvals    | total number of evaluations of the
 *                            | system function F(u) (number of
 *                            | calls made to the user-supplied
 *                            | func routine by the KINSPGMR module
 *                            | member subroutines)
 *                            |
 * KINSpgmrGetLastFlag        | returns the last flag returned by
 *                            | the linear solver
 * -----------------------------------------------------------------
 */

int KINSpgmrGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG);
int KINSpgmrGetNumPrecEvals(void *kinmem, long int *npevals);
int KINSpgmrGetNumPrecSolves(void *kinmem, long int *npsolves);
int KINSpgmrGetNumLinIters(void *kinmem, long int *nliters);
int KINSpgmrGetNumConvFails(void *kinmem, long int *nlcfails);
int KINSpgmrGetNumJtimesEvals(void *kinmem, long int *njvevals);
int KINSpgmrGetNumFuncEvals(void *kinmem, long int *nfevalsSG); 
int KINSpgmrGetLastFlag(void *kinmem, int *flag);

/*
 * -----------------------------------------------------------------
 * KINSpgmrGet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpgmrGet* subroutines
 * are the following:
 *
 * KINSPGMR_SUCCESS : means the routine exited normally [0]
 *
 * KINSPGMR_ILL_INPUT : means at least one input parameter was
 *                      invalid (check error message(s)) [-3]
 *
 * KINSPGMR_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                     given [-1]
 *
 * KINSPGMR_LMEM_NULL : means a NULL KINSPGMR memory block pointer
 *                      was given [-2]
 * -----------------------------------------------------------------
 */

#define KINSPGMR_SUCCESS 0

#define KINSPGMR_MEM_NULL  -1
#define KINSPGMR_LMEM_NULL -2
#define KINSPGMR_ILL_INPUT -3
#define KINSPGMR_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
