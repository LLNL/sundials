/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Header file for the Additive Runge Kutta time step module for ARKode.
 *--------------------------------------------------------------*/

#ifndef _ERKSTEP_H
#define _ERKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <arkode/arkode.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ERKSTEP Constants
  ===============================================================*/



/*===============================================================
  ERKSTEP Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepCreate

  This creates the ERK time step module and attaches it to the 
  main ARKode infrastructure.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepCreate(void* arkode_mem, ARKRhsFn f,
                                  realtype t0, N_Vector y0);

  
/*---------------------------------------------------------------
  ERKStepReInit

  This re-initializes the ERK time step module.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f,
                                  realtype t0, N_Vector y0);

  
/*---------------------------------------------------------------
  ERKStep optional input specification functions -- ALL of these 
  must be called AFTER ERKStepCreate.
-----------------------------------------------------------------
 The following functions can be called to set optional inputs
 to values other than the defaults given below.  

 Function                 |  Optional input / [ default value ]
-----------------------------------------------------------------
 ERKStepSetDefaults       | resets all ERKStep optional inputs to
                          | their default values.  Does not change 
                          | problem-defining function pointers or 
                          | user_data pointer.  Also leaves alone 
                          | any data structures/options related 
                          | to the main ARKode infrastructure 
                          | itself (e.g. root-finding).
                          | [internal]
                          |
 ERKStepSetOrder          | method order to be used by the solver.
                          | [4]
                          |
 ERKStepSetERKTable       | specifies to use a customized Butcher 
                          | table for the explicit portion of the 
                          | system.  This automatically calls 
                          ! ERKStepSetExplicit
                          | [determined by ARKode based on order]
                          |
 ERKStepSetERKTableNum    | specifies to use a built-in Butcher 
                          | table for the explicit portion of the 
                          | system.  The integer argument should 
                          | match an existing method in 
                          | ARKodeLoadButcherTable() within the file
                          | arkode_butcher.c.  Error-checking is
                          | performed to ensure that the table
                          | exists, and is not implicit.  This 
                          ! automatically calls ERKStepSetExplicit
                          | [determined by ARKode based on order]
                          |
 ERKStepSetCFLFraction    | safety factor to use for explicitly 
                          ! stable steps
                          | [0.5]
                          |
 ERKStepSetSafetyFactor   | safety factor to use for error-based 
                          ! step adaptivity
                          | [0.96]
                          |
 ERKStepSetErrorBias      | error bias factor to use in error-based
                          ! step adaptivity
                          | [1.5]
                          |
 ERKStepSetMaxGrowth      | maximum growth factor for successive 
                          ! time steps (not including the first step).
                          | [20.0]
                          |
 ERKStepSetMaxFirstGrowth | maximum growth factor for first step.
                          | [10000.0]
                          |
 ERKStepSetMaxEFailGrowth | maximum growth factor after an error failure.
                          | [0.3]
                          |
 ERKStepSetSmallNumEFails | maximum number of error failures before 
                          ! MaxFailGrowth factor is used.
                          | [2]
                          |
 ERKStepSetFixedStepBounds| step growth interval to force retention of 
                          ! the same step size
                          | [1.0 1.5]
                          |
 ERKStepSetAdaptivityMethod| Method to use for time step adaptivity
                          | [0]
                          |
 ERKStepSetAdaptivityFn   | user-provided time step adaptivity 
                          | function.
                          | [internal]
                          |
 ERKStepSetStabilityFn    | user-provided explicit time step 
                          | stability function.
                          | [internal]
                          |
 ERKStepSetMaxErrTestFails| Maximum number of error test failures
                          | in attempting one step.
                          | [7]
-----------------------------------------------------------------
 Return flag:
   ARK_SUCCESS   if successful
   ARK_MEM_NULL  if the arkode memory is NULL
   ARK_ILL_INPUT if an argument has an illegal value
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ERKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ERKStepSetERKTable(void *arkode_mem, int s, 
                                       int q, int p, realtype *c, 
                                       realtype *A, realtype *b, 
                                       realtype *bembed);
SUNDIALS_EXPORT int ERKStepSetERKTableNum(void *arkode_mem, int itable);
SUNDIALS_EXPORT int ERKStepSetCFLFraction(void *arkode_mem, 
                                          realtype cfl_frac);
SUNDIALS_EXPORT int ERKStepSetSafetyFactor(void *arkode_mem, 
                                           realtype safety);
SUNDIALS_EXPORT int ERKStepSetErrorBias(void *arkode_mem, 
                                        realtype bias);
SUNDIALS_EXPORT int ERKStepSetMaxGrowth(void *arkode_mem, 
                                        realtype mx_growth);
SUNDIALS_EXPORT int ERKStepSetFixedStepBounds(void *arkode_mem, 
                                              realtype lb, realtype ub);
SUNDIALS_EXPORT int ERKStepSetAdaptivityMethod(void *arkode_mem, 
                                               int imethod, 
                                               int idefault, int pq, 
                                               realtype *adapt_params);
SUNDIALS_EXPORT int ERKStepSetAdaptivityFn(void *arkode_mem, 
                                           ARKAdaptFn hfun, 
                                           void *h_data);
SUNDIALS_EXPORT int ERKStepSetMaxFirstGrowth(void *arkode_mem, 
                                             realtype etamx1);
SUNDIALS_EXPORT int ERKStepSetMaxEFailGrowth(void *arkode_mem, 
                                             realtype etamxf);
SUNDIALS_EXPORT int ERKStepSetSmallNumEFails(void *arkode_mem, 
                                             int small_nef);
SUNDIALS_EXPORT int ERKStepSetStabilityFn(void *arkode_mem, 
                                          ARKExpStabFn EStab, 
                                          void *estab_data);
SUNDIALS_EXPORT int ERKStepSetMaxErrTestFails(void *arkode_mem, 
                                              int maxnef);

  
/*---------------------------------------------------------------
  Optional outputs from the ERKStep module:  the following 
  functions can be called to get optional outputs and statistics 
  related to the main integrator.

  ERKStepGetNumExpSteps returns the cumulative number of stability 
                        limited steps taken by the solver

  ERKStepGetNumAccSteps returns the cumulative number of accuracy 
                        limited steps taken by the solver

  ERKStepGetNumStepAttempts returns the total number of steps
                            attempted by the solver

  ERKStepGetNumRhsEvals returns the number of calls to the user's
                        f functions

  ERKStepGetNumErrTestFails returns the number of local error test
                            failures that have occured

  ERKStepGetCurrentButcherTable returns the Butcher table 
                                currently in use

  ERKStepGetEstLocalErrors returns the vector of estimated local
                           errors. The user must allocate space
                           for ele.

  The return value of ERKStepGet* is one of:
     ARK_SUCCESS   if successful
     ARK_MEM_NULL  if the ARKode or ERKStep memory structures 
                   were NULL
     ARK_LMEM_NULL if a linear solver memory structure was NULL
  ---------------------------------------------------------------*/

SUNDIALS_EXPORT int ERKStepGetNumExpSteps(void *arkode_mem, 
                                          long int *expsteps);
SUNDIALS_EXPORT int ERKStepGetNumAccSteps(void *arkode_mem, 
                                          long int *accsteps);
SUNDIALS_EXPORT int ERKStepGetNumStepAttempts(void *arkode_mem, 
                                              long int *step_attempts);
SUNDIALS_EXPORT int ERKStepGetNumRhsEvals(void *arkode_mem, 
                                          long int *nfevals);
SUNDIALS_EXPORT int ERKStepGetNumErrTestFails(void *arkode_mem, 
                                              long int *netfails);
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                                  ARKodeButcherTable *B);
SUNDIALS_EXPORT int ERKStepGetEstLocalErrors(void *arkode_mem, 
                                             N_Vector ele);

/*---------------------------------------------------------------
 As a convenience, the following functions provides the
 optional outputs in one group.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(void *arkode_mem, 
                                               long int *expsteps, 
                                               long int *accsteps, 
                                               long int *step_attempts, 
                                               long int *nfevals, 
                                               long int *netfails);

/*---------------------------------------------------------------
 Function : ERKStepWriteParameters
-----------------------------------------------------------------
 ERKStepWriteParameters outputs all timestepper module parameters
 to the provided file pointer.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
 Function : ERKStepWriteButcher
-----------------------------------------------------------------
 ERKStepWriteButcher outputs the Butcher tables to the 
 provided file pointer.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepWriteButcher(void *arkode_mem, FILE *fp);

  
#ifdef __cplusplus
}
#endif

#endif
