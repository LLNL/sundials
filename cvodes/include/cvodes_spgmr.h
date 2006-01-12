/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-12 22:53:34 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the CVODES scaled preconditioned 
 * GMRES linear solver, CVSPGMR.
 *
 *
 * Part I contains type definitions and function prototypes for using
 * CVSPGMR on forward problems (IVP integration and/or FSA)
 *
 * Part II contains type definitions and function prototypes for using
 * CVSPGMR on adjoint (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _CVsSPGMR_H
#define _CVSSPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_spils.h"
#include "sundials_spgmr.h"

  /*
   * -----------------------------------------------------------------
   * CVSPGMR solver constants
   * -----------------------------------------------------------------
   * CVSPGMR_MAXL   : default value for the maximum Krylov
   *                  dimension
   *
   * CVSPGMR_MSBPRE : maximum number of steps between
   *                  preconditioner evaluations
   *
   * CVSPGMR_DGMAX  : maximum change in gamma between
   *                  preconditioner evaluations
   *
   * CVSPGMR_DELT   : default value for factor by which the
   *                  tolerance on the nonlinear iteration is
   *                  multiplied to get a tolerance on the linear
   *                  iteration
   * -----------------------------------------------------------------
   */

#define CVSPGMR_MAXL   5
#define CVSPGMR_MSBPRE 50
#define CVSPGMR_DGMAX  RCONST(0.2)
#define CVSPGMR_DELT   RCONST(0.05)

  /* CVSPGMR return values */

#define CVSPGMR_SUCCESS    0
#define CVSPGMR_MEM_NULL  -1
#define CVSPGMR_LMEM_NULL -2
#define CVSPGMR_ILL_INPUT -3
#define CVSPGMR_MEM_FAIL  -4

#define CVSPGMR_ADJMEM_NULL -101
#define CVSPGMR_LMEMB_NULL  -102

  /* 
   * -----------------------------------------------------------------
   * PART I - forward problems
   * -----------------------------------------------------------------
   */

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
   *           NONE, LEFT, RIGHT, or BOTH defined in iterative.h.
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
   *    CVSPGMR_SUCCESS   if successful
   *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
   *    CVSPGMR_MEM_FAIL  if there was a memory allocation failure
   *    CVSPGMR_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CVSpgmr(void *cvode_mem, int pretype, int maxl);

  /*
   * -----------------------------------------------------------------
   * Function: CVSpgmrSetPrecType
   * -----------------------------------------------------------------
   * CVSpgmrSetPrecType resets the type of preconditioner, pretype,
   *     from the value set in a prior call to CVSpgmr.
   *     This must be one of NONE, LEFT, RIGHT, or BOTH.
   * -----------------------------------------------------------------
   */

  int CVSpgmrSetPrecType(void *cvode_mem, int pretype);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the CVSPGMR linear solver
   * -----------------------------------------------------------------
   * CVSpgmrSetGSType specifies the type of Gram-Schmidt
   *                orthogonalization to be used. This must be one of
   *                the two enumeration constants MODIFIED_GS or
   *                CLASSICAL_GS defined in iterative.h. These correspond
   *                to using modified Gram-Schmidt and classical
   *                Gram-Schmidt, respectively.
   *                Default value is MODIFIED_GS.
   * CVSpgmrSetDelt specifies the factor by which the tolerance on
   *                the nonlinear iteration is multiplied to get a
   *                tolerance on the linear iteration. This is an
   *                optional input to the CVSPGMR solver.
   *                Default value is 0.05.
   * CVSpgmrSetPreconditioner specifies the PrecSetup and PrecSolve functions.
   *                as well as a pointer to user preconditioner data.
   *                This pointer is passed to PrecSetup and PrecSolve
   *                every time these routines are called.
   *                Default is NULL for al three arguments.
   * CVSpgmrSetJacTimesVecFn specifies the jtimes function and a pointer to
   *                user Jacobian data. This pointer is passed to jtimes every 
   *                time the jtimes routine is called.
   *                Default is to use an internal finite difference
   *                approximation routine.
   *
   * The return value of CVSpgmrSet* is one of:
   *    CVSPGMR_SUCCESS   if successful
   *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
   *    CVSPGMR_LMEM_NULL if the cvspgmr memory was NULL
   *    CVSPGMR_ILL_INPUT if an input has an illegal value
   * -----------------------------------------------------------------
   */

  int CVSpgmrSetGSType(void *cvode_mem, int gstype);
  int CVSpgmrSetDelt(void *cvode_mem, realtype delt);
  int CVSpgmrSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset, 
                               CVSpilsPrecSolveFn psolve, void *P_data);
  int CVSpgmrSetJacTimesVecFn(void *cvode_mem, 
                              CVSpilsJacTimesVecFn jtimes, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the CVSPGMR linear solver
   * -----------------------------------------------------------------
   * CVSpgmrGetWorkSpace returns the real and integer workspace used
   *                     by CVSPGMR.
   * CVSpgmrGetNumPrecEvals returns the number of preconditioner
   *                        evaluations, i.e. the number of calls made
   *                        to PrecSetup with jok==FALSE.
   * CVSpgmrGetNumPrecSolves returns the number of calls made to
   *                         PrecSolve.
   * CVSpgmrGetNumLinIters returns the number of linear iterations.
   * CVSpgmrGetNumConvFails returns the number of linear
   *                        convergence failures.
   * CVSpgmrGetNumJtimesEvals returns the number of calls to jtimes.
   * CVSpgmrGetNumRhsEvals returns the number of calls to the user
   *                       f routine due to finite difference Jacobian
   *                       times vector evaluation.
   * CVSpgmrGetLastFlag returns the last error flag set by any of
   *                    the CVSPGMR interface functions.
   *
   * The return value of CVSpgmrGet* is one of:
   *    CVSPGMR_SUCCESS   if successful
   *    CVSPGMR_MEM_NULL  if the cvode memory was NULL
   *    CVSPGMR_LMEM_NULL if the cvspgmr memory was NULL
   * -----------------------------------------------------------------
   */

  int CVSpgmrGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
  int CVSpgmrGetNumPrecEvals(void *cvode_mem, long int *npevals);
  int CVSpgmrGetNumPrecSolves(void *cvode_mem, long int *npsolves);
  int CVSpgmrGetNumLinIters(void *cvode_mem, long int *nliters);
  int CVSpgmrGetNumConvFails(void *cvode_mem, long int *nlcfails);
  int CVSpgmrGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
  int CVSpgmrGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS); 
  int CVSpgmrGetLastFlag(void *cvode_mem, int *flag);


  /* 
   * -----------------------------------------------------------------
   * PART II - backward problems
   * -----------------------------------------------------------------
   */

  int CVSpgmrB(void *cvadj_mem, int pretypeB, int maxlB);
  
  int CVSpgmrSetPrecTypeB(void *cvadj_mem, int pretypeB);
  int CVSpgmrSetGSTypeB(void *cvadj_mem, int gstypeB);
  int CVSpgmrSetDeltB(void *cvadj_mem, realtype deltB);
  int CVSpgmrSetPreconditionerB(void *cvadj_mem, CVSpilsPrecSetupFnB psetB,
                                CVSpilsPrecSolveFnB psolveB, void *P_dataB);
  int CVSpgmrSetJacTimesVecFnB(void *cvadj_mem, CVSpilsJacTimesVecFnB jtimesB,
                               void *jac_dataB);


#ifdef __cplusplus
}
#endif

#endif
