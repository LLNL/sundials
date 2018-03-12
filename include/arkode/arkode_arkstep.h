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
 * Header file for the Additive Runge Kutta time step module for ARKODE.
 *--------------------------------------------------------------*/

#ifndef _ARKSTEP_H
#define _ARKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKSTEP Constants
  ===============================================================*/



/*===============================================================
  ARKSTEP Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepCreate

  This creates the ARK time step module and attaches it to the 
  main ARKode infrastructure.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepCreate(void* arkode_mem, ARKRhsFn fe,
                                  ARKRhsFn fi, realtype t0, N_Vector y0);

  
/*---------------------------------------------------------------
  ARKStepReInit

  This re-initializes the ARK time step module.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                                  ARKRhsFn fi, realtype t0, N_Vector y0);

  
/*---------------------------------------------------------------
  ARKStep optional input specification functions -- ALL of these 
  must be called AFTER ARKStepCreate.
-----------------------------------------------------------------
 The following functions can be called to set optional inputs
 to values other than the defaults given below.  

 Function                 |  Optional input / [ default value ]
-----------------------------------------------------------------
 ARKStepSetDefaults       | resets all ARKStep optional inputs to
                          | their default values.  Does not change 
                          | problem-defining function pointers or 
                          | user_data pointer.  Also leaves alone 
                          | any data structures/options related 
                          | to the main ARKode infrastructure 
                          | itself (e.g. root-finding).
                          | [internal]
                          |
 ARKStepSetOptimalParams  | sets all adaptivity and solver 
                          | parameters to our 'best guess' values, 
                          | for a given integration method (ERK, 
                          | DIRK, ARK) and a given method order.  
                          | Should only be called after the method
                          | order and integration method have been 
                          ! set.
                          | [internal]
                          |
 ARKStepSetOrder          | method order to be used by the solver.
                          | [4]
                          |
 ARKStepSetLinear         | specifies that the implicit portion of 
                          | the problem is linear, and to tighten 
                          | the linear solver tolerances while 
                          | taking only one Newton iteration.
                          | [SUNFALSE]
                          |
 ARKStepSetNonlinear      | specifies that the implicit portion of 
                          | the problem is nonlinear.  Used to undo
                          | a previous call to ARKStepSetLinear
                          | [SUNTRUE]
                          |
 ARKStepSetFixedPoint     | specifies that the implicit portion of 
                          | the problem should use the accelerated 
                          | fixed-point solver.
                          | [SUNFALSE]
                          |
 ARKStepSetNewton         | specifies that the implicit portion of 
                          | the problem should use the modified Newton 
                          | solver.  Used to undo a previous call to
                          | ARKStepSetFixedPoint
                          | [SUNTRUE]
                          |
 ARKStepSetExplicit       | specifies that implicit portion of 
                          | problem is disabled, and to use an 
                          | explicit RK method.
                          | [SUNFALSE]
                          |
 ARKStepSetImplicit       | specifies that explicit portion of 
                          | problem is disabled, and to use an 
                          | implicit RK method.
                          | [SUNFALSE]
                          |
 ARKStepSetImEx           | specifies that problem has both 
                          | implicit and explicit parts, and to 
                          | use an ARK method.
                          | [SUNTRUE]
                          |
 ARKStepSetERKTable       | specifies to use a customized Butcher 
                          | table for the explicit portion of the 
                          | system.  This automatically calls 
                          ! ARKStepSetExplicit
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetIRKTable       | specifies to use a customized Butcher 
                          | table for the implicit portion of the 
                          | system. This automatically calls 
                          ! ARKStepSetImplicit
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetARKTables      | specifies to use customized Butcher 
                          | tables for the IMEX system.  This 
                          ! automatically calls ARKStepSetImEx
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetERKTableNum    | specifies to use a built-in Butcher 
                          | table for the explicit portion of the 
                          | system.  The integer argument should 
                          | match an existing method in 
                          | ARKodeLoadButcherTable() within the file
                          | arkode_butcher.c.  Error-checking is
                          | performed to ensure that the table
                          | exists, and is not implicit.  This 
                          ! automatically calls ARKStepSetExplicit
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetIRKTableNum    | specifies to use a built-in Butcher 
                          | table for the implicit portion of the 
                          | system.  The integer argument should 
                          | match an existing method in 
                          | ARKodeLoadButcherTable() within the file
                          | arkode_butcher.c.  Error-checking is
                          | performed to ensure that the table
                          | exists, and is not explicit.  This 
                          ! automatically calls ARKStepSetImplicit
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetARKTableNum    | specifies to use a built-in Butcher 
                          | tables for the ImEx system.  The 
                          ! integer arguments should match existing 
                          | methods in ARKodeLoadButcherTable() 
                          | within the file arkode_butcher.c.  
                          | Error-checking is performed to ensure 
                          | that the tables exist.  Subsequent 
                          | error-checking is automatically performed
                          | to ensure that the tables' stage times 
                          | and solution coefficients match.  This 
                          ! automatically calls ARKStepSetImEx
                          | [determined by ARKODE based on order]
                          |
 ARKStepSetCFLFraction    | safety factor to use for explicitly 
                          ! stable steps
                          | [0.5]
                          |
 ARKStepSetSafetyFactor   | safety factor to use for error-based 
                          ! step adaptivity
                          | [0.96]
                          |
 ARKStepSetErrorBias      | error bias factor to use in error-based
                          ! step adaptivity
                          | [1.5]
                          |
 ARKStepSetMaxGrowth      | maximum growth factor for successive 
                          ! time steps (not including the first step).
                          | [20.0]
                          |
 ARKStepSetMaxFirstGrowth | maximum growth factor for first step.
                          | [10000.0]
                          |
 ARKStepSetMaxEFailGrowth | maximum growth factor after an error failure.
                          | [0.3]
                          |
 ARKStepSetSmallNumEFails | maximum number of error failures before 
                          ! MaxFailGrowth factor is used.
                          | [2]
                          |
 ARKStepSetMaxCFailGrowth | maximum growth factor after a convergence failure.
                          | [0.25]
                          |
 ARKStepSetFixedStepBounds| step growth interval to force retention of 
                          ! the same step size
                          | [1.0 1.5]
                          |
 ARKStepSetAdaptivityMethod| Method to use for time step adaptivity
                          | [0]
                          |
 ARKStepSetAdaptivityFn   | user-provided time step adaptivity 
                          | function.
                          | [internal]
                          |
 ARKStepSetNonlinCRDown   | user-provided nonlinear convergence
                          | rate constant.
                          | [0.3]
                          |
 ARKStepSetNonlinRDiv     | user-provided nonlinear divergence ratio.
                          | [2.3]
                          |
 ARKStepSetDeltaGammaMax  | user-provided linear setup decision
                          | constant.
                          | [0.2]
                          |
 ARKStepSetMaxStepsBetweenLSet| user-provided linear setup decision
                          | constant.
                          | [20]
                          |
 ARKStepSetPredictorMethod| Method to use for predicting implicit 
                          | solutions.
                          | [0]
                          |
 ARKStepSetStabilityFn    | user-provided explicit time step 
                          | stability function.
                          | [internal]
                          |
 ARKStepSetMaxErrTestFails| Maximum number of error test failures
                          | in attempting one step.
                          | [7]
                          |
 ARKStepSetMaxNonlinIters | Maximum number of nonlinear solver
                          | iterations at one stage solution.
                          | [3]
                          |
 ARKStepSetMaxConvFails   | Maximum number of convergence failures
                          | allowed in attempting one step.
                          | [10]
                          |
 ARKStepSetNonlinConvCoef | Coefficient in the nonlinear
                          | convergence test.
                          | [0.1]
-----------------------------------------------------------------
 Return flag:
   ARK_SUCCESS   if successful
   ARK_MEM_NULL  if the arkode memory is NULL
   ARK_ILL_INPUT if an argument has an illegal value
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ARKStepSetOptimalParams(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ARKStepSetLinear(void *arkode_mem, int timedepend);
SUNDIALS_EXPORT int ARKStepSetNonlinear(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetFixedPoint(void *arkode_mem, long int fp_m);
SUNDIALS_EXPORT int ARKStepSetNewton(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetExplicit(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImplicit(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImEx(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetERKTable(void *arkode_mem, int s, 
                                       int q, int p, realtype *c, 
                                       realtype *A, realtype *b, 
                                       realtype *bembed);
SUNDIALS_EXPORT int ARKStepSetIRKTable(void *arkode_mem, int s, 
                                       int q, int p, realtype *c, 
                                       realtype *A, realtype *b, 
                                       realtype *bembed);
SUNDIALS_EXPORT int ARKStepSetARKTables(void *arkode_mem, int s, 
                                        int q, int p, 
                                        realtype *ci, realtype *ce, 
                                        realtype *Ai, realtype *Ae, 
                                        realtype *bi, realtype *be, 
                                        realtype *b2i, realtype *b2e);
SUNDIALS_EXPORT int ARKStepSetERKTableNum(void *arkode_mem, int itable);
SUNDIALS_EXPORT int ARKStepSetIRKTableNum(void *arkode_mem, int itable);
SUNDIALS_EXPORT int ARKStepSetARKTableNum(void *arkode_mem, 
                                          int itable, int etable);
SUNDIALS_EXPORT int ARKStepSetCFLFraction(void *arkode_mem, 
                                          realtype cfl_frac);
SUNDIALS_EXPORT int ARKStepSetSafetyFactor(void *arkode_mem, 
                                           realtype safety);
SUNDIALS_EXPORT int ARKStepSetErrorBias(void *arkode_mem, 
                                        realtype bias);
SUNDIALS_EXPORT int ARKStepSetMaxGrowth(void *arkode_mem, 
                                        realtype mx_growth);
SUNDIALS_EXPORT int ARKStepSetFixedStepBounds(void *arkode_mem, 
                                              realtype lb, realtype ub);
SUNDIALS_EXPORT int ARKStepSetAdaptivityMethod(void *arkode_mem, 
                                               int imethod, 
                                               int idefault, int pq, 
                                               realtype *adapt_params);
SUNDIALS_EXPORT int ARKStepSetAdaptivityFn(void *arkode_mem, 
                                           ARKAdaptFn hfun, 
                                           void *h_data);
SUNDIALS_EXPORT int ARKStepSetMaxFirstGrowth(void *arkode_mem, 
                                             realtype etamx1);
SUNDIALS_EXPORT int ARKStepSetMaxEFailGrowth(void *arkode_mem, 
                                             realtype etamxf);
SUNDIALS_EXPORT int ARKStepSetSmallNumEFails(void *arkode_mem, 
                                             int small_nef);
SUNDIALS_EXPORT int ARKStepSetMaxCFailGrowth(void *arkode_mem, 
                                             realtype etacf);
SUNDIALS_EXPORT int ARKStepSetNonlinCRDown(void *arkode_mem, 
                                           realtype crdown);
SUNDIALS_EXPORT int ARKStepSetNonlinRDiv(void *arkode_mem, 
                                         realtype rdiv);
SUNDIALS_EXPORT int ARKStepSetDeltaGammaMax(void *arkode_mem, 
                                            realtype dgmax);
SUNDIALS_EXPORT int ARKStepSetMaxStepsBetweenLSet(void *arkode_mem, 
                                                  int msbp);
SUNDIALS_EXPORT int ARKStepSetPredictorMethod(void *arkode_mem, 
                                              int method);
SUNDIALS_EXPORT int ARKStepSetStabilityFn(void *arkode_mem, 
                                          ARKExpStabFn EStab, 
                                          void *estab_data);
SUNDIALS_EXPORT int ARKStepSetMaxErrTestFails(void *arkode_mem, 
                                              int maxnef);
SUNDIALS_EXPORT int ARKStepSetMaxNonlinIters(void *arkode_mem, 
                                             int maxcor);
SUNDIALS_EXPORT int ARKStepSetMaxConvFails(void *arkode_mem, 
                                           int maxncf);
SUNDIALS_EXPORT int ARKStepSetNonlinConvCoef(void *arkode_mem, 
                                             realtype nlscoef);

  
/*---------------------------------------------------------------
  Optional outputs from the ARKStep module:  the following 
  functions can be called to get optional outputs and statistics 
  related to the main integrator.

  ARKStepGetNumExpSteps returns the cumulative number of stability 
                        limited steps taken by the solver

  ARKStepGetNumAccSteps returns the cumulative number of accuracy 
                        limited steps taken by the solver

  ARKStepGetNumStepAttempts returns the total number of steps
                            attempted by the solver

  ARKStepGetNumRhsEvals returns the number of calls to the user's
                        fe and fi functions

  ARKStepGetNumLinSolvSetups returns the number of calls made to
                             the linear solver's setup routine

  ARKStepGetNumMassSolves returns the number of calls made to
                          the mass matrix solve routine

  ARKStepGetNumMassMultiplies returns the number of calls made to
                              the mass matrix times vector routine

  ARKStepGetNumErrTestFails returns the number of local error test
                            failures that have occured

  ARKStepGetCurrentButcherTables returns the explicit and implicit
                                 Butcher tables currently in use

  ARKStepGetEstLocalErrors returns the vector of estimated local
                           errors. The user must allocate space
                           for ele.

  The return value of ARKStepGet* is one of:
     ARK_SUCCESS   if successful
     ARK_MEM_NULL  if the ARKode or ARKStep memory structures 
                   were NULL
     ARK_LMEM_NULL if a linear solver memory structure was NULL
  ---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKStepGetNumExpSteps(void *arkode_mem, 
                                          long int *expsteps);
SUNDIALS_EXPORT int ARKStepGetNumAccSteps(void *arkode_mem, 
                                          long int *accsteps);
SUNDIALS_EXPORT int ARKStepGetNumStepAttempts(void *arkode_mem, 
                                              long int *step_attempts);
SUNDIALS_EXPORT int ARKStepGetNumRhsEvals(void *arkode_mem, 
                                          long int *nfe_evals, 
                                          long int *nfi_evals);
SUNDIALS_EXPORT int ARKStepGetNumLinSolvSetups(void *arkode_mem, 
                                               long int *nlinsetups);
SUNDIALS_EXPORT int ARKStepGetNumMassSolves(void *arkode_mem, 
                                            long int *nMassSolves);
SUNDIALS_EXPORT int ARKStepGetNumMassMultiplies(void *arkode_mem, 
                                                long int *nMassMult);
SUNDIALS_EXPORT int ARKStepGetNumErrTestFails(void *arkode_mem, 
                                              long int *netfails);
SUNDIALS_EXPORT int ARKStepGetCurrentButcherTables(void *arkode_mem,
                                                   ARKodeButcherTable *Bi,
                                                   ARKodeButcherTable *Be);
SUNDIALS_EXPORT int ARKStepGetEstLocalErrors(void *arkode_mem, 
                                             N_Vector ele);

/*---------------------------------------------------------------
 As a convenience, the following functions provides the
 optional outputs in one group.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetTimestepperStats(void *arkode_mem, 
                                               long int *expsteps, 
                                               long int *accsteps, 
                                               long int *step_attempts, 
                                               long int *nfe_evals, 
                                               long int *nfi_evals, 
                                               long int *nlinsetups, 
                                               long int *netfails);

/*---------------------------------------------------------------
 Nonlinear solver optional output extraction functions
-----------------------------------------------------------------
 The following functions can be called to get optional outputs
 and statistics related to the nonlinear solver.
-----------------------------------------------------------------
 ARKStepGetNumNonlinSolvIters returns the number of nonlinear
                             solver iterations performed.

 ARKStepGetNumNonlinSolvConvFails returns the number of nonlinear
                                 convergence failures.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetNumNonlinSolvIters(void *arkode_mem, 
                                                long int *nniters);
SUNDIALS_EXPORT int ARKStepGetNumNonlinSolvConvFails(void *arkode_mem, 
                                                    long int *nncfails);

/*---------------------------------------------------------------
 As a convenience, the following function provides the
 nonlinear solver optional outputs in a group.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetNonlinSolvStats(void *arkode_mem, 
                                             long int *nniters,
                                             long int *nncfails);

/*---------------------------------------------------------------
 Function : ARKStepWriteParameters
-----------------------------------------------------------------
 ARKStepWriteParameters outputs all timestepper module parameters
 to the provided file pointer.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
 Function : ARKStepWriteButcher
-----------------------------------------------------------------
 ARKStepWriteButcher outputs the Butcher tables to the 
 provided file pointer.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepWriteButcher(void *arkode_mem, FILE *fp);

  
#ifdef __cplusplus
}
#endif

#endif
