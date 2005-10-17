/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2005-10-17 21:50:04 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES scaled
 * preconditioned Bi-CGSTAB linear solver, CVSPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPBCG_H
#define _CVSPBCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvspils.h"
#include "spbcg.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * CVSPBCG solver constants
 * -----------------------------------------------------------------
 * CVSPBCG_MAXL   : default value for the maximum Krylov
 *                  dimension
 *
 * CVSPBCG_MSBPRE : maximum number of steps between
 *                  preconditioner evaluations
 *
 * CVSPBCG_DGMAX  : maximum change in gamma between
 *                  preconditioner evaluations
 *
 * CVSPBCG_DELT   : default value for factor by which the
 *                  tolerance on the nonlinear iteration is
 *                  multiplied to get a tolerance on the linear
 *                  iteration
 * -----------------------------------------------------------------
 */

#define CVSPBCG_MAXL   5
#define CVSPBCG_MSBPRE 50
#define CVSPBCG_DGMAX  RCONST(0.2)
#define CVSPBCG_DELT   RCONST(0.05)

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
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_MEM_FAIL  if there was a memory allocation failure
 *    CVSPBCG_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int CVSpbcg(void *cvode_mem, int pretype, int maxl);

/*
 * -----------------------------------------------------------------
 * Function: CVSpbcgSetPrecType
 * -----------------------------------------------------------------
 * CVSpbcgSetPrecType resets the type of preconditioner, pretype,
 *                    from the value set in a prior call to CVSpbcg.
 *                    This must be one of PREC_NONE, PREC_LEFT,
 *                    PREC_RIGHT, or PREC_BOTH.
 * -----------------------------------------------------------------
 */

int CVSpbcgSetPrecType(void *cvode_mem, int pretype);

/*
 * -----------------------------------------------------------------
 * Function: CVSpbcgSetMaxl
 * -----------------------------------------------------------------
 * CVSpbcgSetMaxl     resets the maximum Krylov subspace size, maxl,
 *                    from the value set in a prior call to CVSpbcg.
 *                    An input value <= 0, gives the default value.
 * -----------------------------------------------------------------
 */

int CVSpbcgSetMaxl(void *cvode_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPBCG linear solver
 * -----------------------------------------------------------------
 * CVSpbcgSetDelt specifies the factor by which the tolerance on
 *                   the nonlinear iteration is multiplied to get a
 *                   tolerance on the linear iteration. This is an
 *                   optional input to the CVSPBCG solver.
 *                   Default value is 0.05.
 * CVSpbcgSetPreconditionr specifies the PrecSetup and PrecSolve
 *                   functions and the pointer to user data that
 *                   is pased to these functinos whenever they
 *                   are called.
 *                   Default is NULL for all three.
 * CVSpbcgSetJacTimesVecFn specifies the jtimes function and the 
 *                   pointer to user data that is passed to jtimes
 *                   whenever it is called.
 *                   Default is to use an internal finite
 *                   difference approximation routine.
 *
 * The return value of CVSpbcgSet* is one of:
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_LMEM_NULL if the cvspbcg memory was NULL
 *    CVSPBCG_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

int CVSpbcgSetDelt(void *cvode_mem, realtype delt);
int CVSpbcgSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset,
                             CVSpilsPrecSolveFn psolve, void *P_data);
int CVSpbcgSetJacTimesVecFn(void *cvode_mem, 
                            CVSpilsJacTimesVecFn jtimes, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSPBCG linear solver
 * -----------------------------------------------------------------
 * CVSpbcgGetWorkSpace returns the real and integer work space used
 *                     by CVSPBCG.
 * CVSpbcgGetNumPrecEvals returns the number of preconditioner
 *                        evaluations, i.e. the number of calls made
 *                        to PrecSetup with jok==FALSE.
 * CVSpbcgGetNumPrecSolves returns the number of calls made to
 *                         PrecSolve.
 * CVSpbcgGetNumLinIters returns the number of linear iterations.
 * CVSpbcgGetNumConvFails returns the number of linear
 *                        convergence failures.
 * CVSpbcgGetNumJtimesEvals returns the number of calls to jtimes.
 * CVSpbcgGetNumRhsEvals returns the number of calls to the user
 *                       f routine due to finite difference Jacobian
 *                       times vector evaluation.
 * CVSpbcgGetLastFlag returns the last error flag set by any of
 *                    the CVSPBCG interface functions.
 *
 * The return value of CVSpbcgGet* is one of:
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_LMEM_NULL if the cvspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int CVSpbcgGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
int CVSpbcgGetNumPrecEvals(void *cvode_mem, long int *npevals);
int CVSpbcgGetNumPrecSolves(void *cvode_mem, long int *npsolves);
int CVSpbcgGetNumLinIters(void *cvode_mem, long int *nliters);
int CVSpbcgGetNumConvFails(void *cvode_mem, long int *nlcfails);
int CVSpbcgGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
int CVSpbcgGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS); 
int CVSpbcgGetLastFlag(void *cvode_mem, int *flag);

/* CVSPBCG return values */

#define CVSPBCG_SUCCESS    0
#define CVSPBCG_MEM_NULL  -1
#define CVSPBCG_LMEM_NULL -2
#define CVSPBCG_ILL_INPUT -3
#define CVSPBCG_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
