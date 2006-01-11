/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:56 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * This is the public header file for the IDAS scaled preconditioned
 * TFQMR linear solver module, IDASPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPTFQMR_H
#define _IDASPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idas_spils.h"
#include "sundials_sptfqmr.h"

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
 *    IDASPTFQMR_SUCCESS    if successful
 *    IDASPTFQMR_MEM_NULL   if the ida memory was NULL
 *    IDASPTFQMR_MEM_FAIL   if there was a memory allocation failure
 *    IDASPTFQMR_ILL_INPUT  if there was illegal input.
 * -----------------------------------------------------------------
 */

int IDASptfqmr(void *ida_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDASPTFQMR linear solver
 * -----------------------------------------------------------------
 * IDASptfqmrSetPreconditionr specifies the PrecSetup and PrecSolve
 *                   functions and the pointer to user data that
 *                   is pased to these functinos whenever they
 *                   are called.
 *                   [NULL/NULL/NULL]
 * IDASptfqmrSetJacTimesVecFn specifies the jtimes function and the 
 *                   pointer to user data that is passed to jtimes
 *                   whenever it is called.
 *                   [IDASptfqmrDQJtimes/NULL]
 * IDASptfqmrSetEpsLin specifies the factor in the linear iteration
 *                   convergence test constant.
 *                   [0.05]
 * IDASptfqmrSetIncrementFactor specifies a factor in the increments
 *                   to yy used in the difference quotient
 *                   approximations to matrix-vector products Jv.
 *                   [1.0]
 *
 * The return value of IDASptfqmrSet* is one of:
 *    IDASPTFQMR_SUCCESS   if successful
 *    IDASPTFQMR_MEM_NULL  if the ida memory was NULL
 *    IDASPTFQMR_LMEM_NULL if the idasptfqmr memory was NULL
 * -----------------------------------------------------------------
 */

int IDASptfqmrSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
				IDASpilsPrecSolveFn psolve, void *prec_data);
int IDASptfqmrSetJacTimesVecFn(void *ida_mem, 
			       IDASpilsJacTimesVecFn jtimes, void *jac_data);
int IDASptfqmrSetEpsLin(void *ida_mem, realtype eplifac);
int IDASptfqmrSetIncrementFactor(void *ida_mem, realtype dqincfac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDASPTFQMR linear solver
 * -----------------------------------------------------------------
 * IDASptfqmrGetWorkSpace returns the real and integer workspace used
 *                        by IDASPTFQMR.
 * IDASptfqmrGetNumPrecEvals returns the number of preconditioner
 *                        evaluations, i.e., the number of calls made
 *                        to PrecSetup with jok==FALSE.
 * IDASptfqmrGetNumPrecSolves returns the number of calls made to
 *                        PrecSolve.
 * IDASptfqmrGetNumLinIters returns the number of linear iterations.
 * IDASptfqmrGetNumConvFails returns the number of linear
 *                        convergence failures.
 * IDASptfqmrGetNumJtimesEvals returns the number of calls to jtimes
 * IDASptfqmrGetNumResEvals returns the number of calls to the user
 *                        res routine due to finite difference Jacobian
 *                        times vector evaluation.
 * IDASptfqmrGetLastFlag returns the last error flag set by any of
 *                        the IDASPTFQMR interface functions.
 *
 * The return value of IDASptfqmrGet* is one of:
 *    IDASPTFQMR_SUCCESS   if successful
 *    IDASPTFQMR_MEM_NULL  if the ida memory was NULL
 *    IDASPTFQMR_LMEM_NULL if the idasptfqmr memory was NULL
 * -----------------------------------------------------------------
 */

int IDASptfqmrGetWorkSpace(void *ida_mem, long int *lenrwSG, long int *leniwSG);
int IDASptfqmrGetNumPrecEvals(void *ida_mem, long int *npevals);
int IDASptfqmrGetNumPrecSolves(void *ida_mem, long int *npsolves);
int IDASptfqmrGetNumLinIters(void *ida_mem, long int *nliters);
int IDASptfqmrGetNumConvFails(void *ida_mem, long int *nlcfails);
int IDASptfqmrGetNumJtimesEvals(void *ida_mem, long int *njvevals);
int IDASptfqmrGetNumResEvals(void *ida_mem, long int *nrevalsSG); 
int IDASptfqmrGetLastFlag(void *ida_mem, int *flag);

/* IDASPTFQMR return values */

#define IDASPTFQMR_SUCCESS    0
#define IDASPTFQMR_MEM_NULL  -1 
#define IDASPTFQMR_LMEM_NULL -2 
#define IDASPTFQMR_ILL_INPUT -3
#define IDASPTFQMR_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
