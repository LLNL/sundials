/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-10-03 22:39:45 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES scaled
 * preconditioned TFQMR linear solver, CVSPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPTFQMR_H
#define _CVSPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvspils.h"
#include "sptfqmr.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * CVSPTFQMR solver constants
 * -----------------------------------------------------------------
 * CVSPTFQMR_MAXL   : default value for the maximum Krylov
 *                    dimension
 *
 * CVSPTFQMR_MSBPRE : maximum number of steps between
 *                    preconditioner evaluations
 *
 * CVSPTFQMR_DGMAX  : maximum change in gamma between
 *                    preconditioner evaluations
 *
 * CVSPTFQMR_DELT   : default value for factor by which the
 *                    tolerance on the nonlinear iteration is
 *                    multiplied to get a tolerance on the linear
 *                    iteration
 * -----------------------------------------------------------------
 */

#define CVSPTFQMR_MAXL   5
#define CVSPTFQMR_MSBPRE 50
#define CVSPTFQMR_DGMAX  RCONST(0.2)
#define CVSPTFQMR_DELT   RCONST(0.05)

/*
 * -----------------------------------------------------------------
 * Function : CVSptfqmr
 * -----------------------------------------------------------------
 * A call to the CVSptfqmr function links the main CVODE integrator
 * with the CVSPTFQMR linear solver.
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
 *           optional input to the CVSPTFQMR solver. Pass 0 to
 *           use the default value CVSPTFQMR_MAXL=5.
 *
 * The return value of CVSptfqmr is one of:
 *    CVSPTFQMR_SUCCESS   if successful
 *    CVSPTFQMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPTFQMR_MEM_FAIL  if there was a memory allocation failure
 *    CVSPTFQMR_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int CVSptfqmr(void *cvode_mem, int pretype, int maxl);

/*
 * -----------------------------------------------------------------
 * Function: CVSptfqmrSetPrecType
 * -----------------------------------------------------------------
 * CVSptfqmrSetPrecType resets the type of preconditioner, pretype,
 *                      from the value set in a prior call to CVSptfqmr.
 *                      This must be one of PREC_NONE, PREC_LEFT,
 *                      PREC_RIGHT, or PREC_BOTH.
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetPrecType(void *cvode_mem, int pretype);

/*
 * -----------------------------------------------------------------
 * Function: CVSptfqmrSetMaxl
 * -----------------------------------------------------------------
 * CVSptfqmrSetMaxl   resets the maximum Krylov subspace size, maxl,
 *                    from the value set in a prior call to CVSptfqmr.
 *                    An input value <= 0, gives the default value.
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetMaxl(void *cvode_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPTFQMR linear solver
 * -----------------------------------------------------------------
 * CVSptfqmrSetDelt specifies the factor by which the tolerance on
 *                  the nonlinear iteration is multiplied to get a
 *                  tolerance on the linear iteration. This is an
 *                  optional input to the CVSPTFQMR solver.
 *                  [0.05]
 * CVSptfqmrSetPreconditionr specifies the PrecSetup and PrecSolve
 *                  functions and the pointer to user data that
 *                  is pased to these functinos whenever they
 *                  are called.
 *                  [NULL/NULL/NULL]
 * CVSptfqmrSetJacTimesVecFn specifies the jtimes function and the 
 *                  pointer to user data that is passed to jtimes
 *                  whenever it is called.
 *                  [CVSptfqmrDQJtimes/NULL]
 *
 * The return value of CVSptfqmrSet* is one of:
 *    CVSPTFQMR_SUCCESS   if successful
 *    CVSPTFQMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPTFQMR_LMEM_NULL if the cvsptfqmr memory was NULL
 *    CVSPTFQMR_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

int CVSptfqmrSetDelt(void *cvode_mem, realtype delt);
int CVSptfqmrSetPreconditioner(void *cvode_mem,
			       CVSpilsPrecSetupFn pset,
			       CVSpilsPrecSolveFn psolve,
			       void *P_data);
int CVSptfqmrSetJacTimesVecFn(void *cvode_mem, 
			      CVSpilsJacTimesVecFn jtimes,
			      void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSPTFQMR linear solver
 * -----------------------------------------------------------------
 * CVSptfqmrGetWorkSpace returns the real and integer work space used
 *                       by CVSPTFQMR.
 * CVSptfqmrGetNumPrecEvals returns the number of preconditioner
 *                          evaluations, i.e., the number of calls made
 *                          to PrecSetup with jok==FALSE.
 * CVSptfqmrGetNumPrecSolves returns the number of calls made to
 *                           PrecSolve.
 * CVSptfqmrGetNumLinIters returns the number of linear iterations.
 * CVSptfqmrGetNumConvFails returns the number of linear
 *                          convergence failures.
 * CVSptfqmrGetNumJtimesEvals returns the number of calls to jtimes.
 * CVSptfqmrGetNumRhsEvals returns the number of calls to the user
 *                         f routine due to finite difference Jacobian
 *                         times vector evaluation.
 * CVSptfqmrGetLastFlag returns the last error flag set by any of
 *                      the CVSPTFQMR interface functions.
 *
 * The return value of CVSptfqmrGet* is one of:
 *    CVSPTFQMR_SUCCESS   if successful
 *    CVSPTFQMR_MEM_NULL  if the cvode memory was NULL
 *    CVSPTFQMR_LMEM_NULL if the cvsptfqmr memory was NULL
 * -----------------------------------------------------------------
 */

int CVSptfqmrGetWorkSpace(void *cvode_mem, long int *lenrwSG, long int *leniwSG);
int CVSptfqmrGetNumPrecEvals(void *cvode_mem, long int *npevals);
int CVSptfqmrGetNumPrecSolves(void *cvode_mem, long int *npsolves);
int CVSptfqmrGetNumLinIters(void *cvode_mem, long int *nliters);
int CVSptfqmrGetNumConvFails(void *cvode_mem, long int *nlcfails);
int CVSptfqmrGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
int CVSptfqmrGetNumRhsEvals(void *cvode_mem, long int *nfevalsSG); 
int CVSptfqmrGetLastFlag(void *cvode_mem, int *flag);

/* CVSPTFQMR return values */

#define CVSPTFQMR_SUCCESS    0
#define CVSPTFQMR_MEM_NULL  -1
#define CVSPTFQMR_LMEM_NULL -2
#define CVSPTFQMR_ILL_INPUT -3
#define CVSPTFQMR_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
