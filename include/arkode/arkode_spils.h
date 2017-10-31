/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * Header file for the Scaled, Preconditioned Iterative Linear 
 * Solver interface in ARKODE.
 *--------------------------------------------------------------*/

#ifndef _ARKSPILS_H
#define _ARKSPILS_H

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
 ARKSPILS Constants
===============================================================*/

/*---------------------------------------------------------------
 ARKSPILS return values -- ADJUST CONSTANTS AS NECESSARY
---------------------------------------------------------------*/
#define ARKSPILS_SUCCESS       0
#define ARKSPILS_MEM_NULL     -1
#define ARKSPILS_LMEM_NULL    -2
#define ARKSPILS_ILL_INPUT    -3
#define ARKSPILS_MEM_FAIL     -4
#define ARKSPILS_PMEM_NULL    -5
#define ARKSPILS_MASSMEM_NULL -6

#define ARKSPILS_SUNLS_FAIL  -10

/*---------------------------------------------------------------
 ARKSPILS solver constants:

 ARKSPILS_MAXL   : default value for the maximum Krylov
                   dimension

 ARKSPILS_MSBPRE : maximum number of steps between
                   preconditioner evaluations

 ARKSPILS_DGMAX  : maximum change in gamma between
                   preconditioner evaluations

 ARKSPILS_EPLIN  : default value for factor by which the
                   tolerance on the nonlinear iteration is
                   multiplied to get a tolerance on the linear
                   iteration
---------------------------------------------------------------*/
#define ARKSPILS_MSBPRE 50 /* -- TURN INTO A PARAMETER, MAKE THIS THE DEFAULT */
#define ARKSPILS_DGMAX  RCONST(0.2) /* -- TURN INTO A PARAMETER, MAKE THIS THE DEFAULT */
#define ARKSPILS_EPLIN  RCONST(0.05)


  
/*===============================================================
 ARKSPILS user-supplied function prototypes
===============================================================*/

/*---------------------------------------------------------------
 Type: ARKSpilsPrecSetupFn

 The user-supplied preconditioner setup function PrecSetup and
 the user-supplied preconditioner solve function PrecSolve
 together must define left and right preconditoner matrices
 P1 and P2 (either of which may be trivial), such that the
 product P1*P2 is an approximation to the Newton matrix
 M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
 and gamma is a scalar proportional to the integration step
 size h.  The solution of systems P z = r, with P = P1 or P2,
 is to be carried out by the PrecSolve function, and PrecSetup
 is to do any necessary setup operations.

 The user-supplied preconditioner setup function PrecSetup
 is to evaluate and preprocess any Jacobian-related data
 needed by the preconditioner solve function PrecSolve.
 This might include forming a crude approximate Jacobian,
 and performing an LU factorization on the resulting
 approximation to M.  This function will not be called in
 advance of every call to PrecSolve, but instead will be called
 only as often as necessary to achieve convergence within the
 Inexact Newton iteration.  If the PrecSolve function needs no
 preparation, the PrecSetup function can be NULL.

 For greater efficiency, the PrecSetup function may save
 Jacobian-related data and reuse it, rather than generating it
 from scratch.  In this case, it should use the input flag jok
 to decide whether to recompute the data, and set the output
 flag *jcurPtr accordingly.

 Each call to the PrecSetup function is preceded by a call to
 the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup
 function can use any auxiliary data that is computed and
 saved by the f function and made accessible to PrecSetup.

 A function PrecSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
          namely the predicted value of y(t).

 fy      is the vector f(t,y).

 jok     is an input flag indicating whether Jacobian-related
         data needs to be recomputed, as follows:
           jok == SUNFALSE means recompute Jacobian-related data
                  from scratch.
           jok == SUNTRUE  means that Jacobian data, if saved from
                  the previous PrecSetup call, can be reused
                  (with the current value of gamma).
         A Precset call with jok == SUNTRUE can only occur after
         a call with jok == SUNFALSE.

 jcurPtr is a pointer to an output integer flag which is
         to be set by PrecSetup as follows:
         Set *jcurPtr = SUNTRUE if Jacobian data was recomputed.
         Set *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
                        but saved data was reused.

 gamma   is the scalar appearing in the Newton matrix.

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

 NOTE: If the user's preconditioner needs other quantities,
       they are accessible as follows: hcur (the current stepsize)
       and ewt (the error weight vector) are accessible through
       ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively).
       The unit roundoff is available as UNIT_ROUNDOFF defined in
       sundials_types.h.

 Returned value:
 The value to be returned by the PrecSetup function is a flag
 indicating whether it was successful.  This value should be
   0   if successful,
   > 0 for a recoverable error (step will be retried),
   < 0 for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsPrecSetupFn)(realtype t, N_Vector y, 
                                   N_Vector fy, booleantype jok, 
                                   booleantype *jcurPtr,
                                   realtype gamma, void *user_data);


/*---------------------------------------------------------------
 Type: ARKSpilsPrecSolveFn

 The user-supplied preconditioner solve function PrecSolve
 is to solve a linear system P z = r in which the matrix P is
 one of the preconditioner matrices P1 or P2, depending on the
 type of preconditioning chosen.

 A function PrecSolve must have the prototype given below.
 Its parameters are as follows:

 t      is the current value of the independent variable.

 y      is the current value of the dependent variable vector.

 fy     is the vector f(t,y).

 r      is the right-hand side vector of the linear system.

 z      is the output vector computed by PrecSolve.

 gamma  is the scalar appearing in the Newton matrix.

 delta  is an input tolerance for use by PSolve if it uses
        an iterative method in its solution.  In that case,
        the residual vector Res = r - P z of the system
        should be made less than delta in weighted L2 norm,
        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
        Note: the error weight vector ewt can be obtained
        through a call to the routine ARKodeGetErrWeights.

 lr     is an input flag indicating whether PrecSolve is to use
        the left preconditioner P1 or right preconditioner
        P2: lr = 1 means use P1, and lr = 2 means use P2.

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

 Returned value:
 The value to be returned by the PrecSolve function is a flag
 indicating whether it was successful.  This value should be
   0 if successful,
   positive for a recoverable error (step will be retried),
   negative for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsPrecSolveFn)(realtype t, N_Vector y, 
                                   N_Vector fy, N_Vector r, 
                                   N_Vector z, realtype gamma, 
                                   realtype delta, int lr, 
                                   void *user_data);


/*---------------------------------------------------------------
 Type: ARKSpilsJacTimesSetupFn

 The user-supplied Jacobian-times-vector product setup function 
 JacTimesSetup and the user-supplied Jacobian-times-vector 
 product function JTimes together must generate the product
 J*v for v, where J is the Jacobian df/dy, or an approximation 
 to it, and v is a given vector. It should return 0 if 
 successful a positive value for a recoverable error or a
 negative value for an unrecoverable failure.

 Each call to the JacTimesSetup function is preceded by a call 
 to the RhsFn fi with the same (t,y) arguments.  Thus the 
 JacTimesSetup function can use any auxiliary data that is 
 computed and saved by the fi function and made accessible to 
 JacTimesSetup.

 A function JacTimesSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
          namely the predicted value of y(t).

 fy      is the vector f(t,y).

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

 Returned value:
 The value to be returned by the JacTimesSetup function is a flag
 indicating whether it was successful.  This value should be
   0   if successful,
   > 0 for a recoverable error (step will be retried),
   < 0 for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsJacTimesSetupFn)(realtype t, N_Vector y, 
                                       N_Vector fy, void *user_data);


/*---------------------------------------------------------------
 Type: ARKSpilsJacTimesVecFn

 The user-supplied function jtimes is to generate the product
 J*v for given v, where J is the Jacobian df/dy, or an
 approximation to it, and v is a given vector. It should return
 0 if successful a positive value for a recoverable error or 
 a negative value for an unrecoverable failure.

 A function jtimes must have the prototype given below. Its
 parameters are as follows:

   v        is the N_Vector to be multiplied by J.

   Jv       is the output N_Vector containing J*v.

   t        is the current value of the independent variable.

   y        is the current value of the dependent variable
            vector.

   fy       is the vector f(t,y).

   user_data   is a pointer to user data, the same as the user_data
            parameter passed to the ARKodeSetUserData function.

   tmp      is a pointer to memory allocated for an N_Vector
            which can be used by Jtimes for work space.
---------------------------------------------------------------*/
typedef int (*ARKSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, 
                                     realtype t, N_Vector y, 
                                     N_Vector fy, void *user_data, 
                                     N_Vector tmp);


/*---------------------------------------------------------------
 Type: ARKSpilsMassTimesSetupFn

 The user-supplied mass matrix-times-vector product setup 
 function MassTimesSetup and the user-supplied mass 
 matrix-times-vector product function MTimes together must 
 generate the product M*v for v, where M is the mass matrix, or 
 an approximation to it, and v is a given vector. It should 
 return 0 if successful a positive value for a recoverable error 
 or a negative value for an unrecoverable failure.

 A function MassTimesSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 mtimes_data  is a pointer to user data - the same as the 
         parameter passed to the ARKodeSetMassTimesVecFn
         function.

 Returned value:
 The value to be returned by the MassTimesSetup function is a 
 flag indicating whether it was successful.  This value should be
   0   if successful,
   > 0 for a recoverable error (step will be retried),
   < 0 for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsMassTimesSetupFn)(realtype t, void *mtimes_data);


/*---------------------------------------------------------------
 Type: ARKSpilsMassTimesVecFn

 The user-supplied function mtimes is to generate the product
 M*v for given v, where M is the mass matrix, or an 
 approximation to it, and v is a given vector. It should return 
 0 if successful or a negative value for an unrecoverable failure.

 A function mtimes must have the prototype given below. Its
 parameters are as follows:

   v        is the N_Vector to be multiplied by M.

   Mv       is the output N_Vector containing M*v.

   t        is the current value of the independent variable.

   mtimes_data   is a pointer to user data, the same as the 
            parameter passed to the ARKodeSetMassTimesVecFn 
            function.
---------------------------------------------------------------*/
typedef int (*ARKSpilsMassTimesVecFn)(N_Vector v, N_Vector Mv, 
                                      realtype t, void *mtimes_data);


/*---------------------------------------------------------------
 Type: ARKSpilsMassPrecSetupFn

 The user-supplied mass matrix preconditioner setup function 
 MPrecSetup and the user-supplied mass matrix preconditioner solve 
 function PrecSolve together must define left and right 
 preconditoner matrices P1 and P2 (either of which may be 
 trivial), such that the product P1*P2 is an approximation to 
 the mass matrix M.  The solution of systems P z = r, with P = P1 
 or P2, is to be carried out by the PrecSolve function, and 
 MPrecSetup is to do any necessary setup operations.

 The user-supplied preconditioner setup function MPrecSetup
 is to evaluate and preprocess any mass-matrix-related data
 needed by the preconditioner solve function PrecSolve.

 A function MPrecSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

 Returned value:
 The value to be returned by the MPrecSetup function is a flag
 indicating whether it was successful.  This value should be
   0   if successful,
   < 0 for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsMassPrecSetupFn)(realtype t, void *user_data);


/*---------------------------------------------------------------
 Type: ARKSpilsMassPrecSolveFn

 The user-supplied mass matrix preconditioner solve function 
 MPrecSolve is to solve a linear system P z = r in which the 
 matrix P is one of the preconditioner matrices P1 or P2, 
 depending on the type of preconditioning chosen.

 A function MPrecSolve must have the prototype given below.
 Its parameters are as follows:

 t      is the current value of the independent variable.

 r      is the right-hand side vector of the linear system.

 z      is the output vector computed by MPrecSolve.

 delta  is an input tolerance for use by PSolve if it uses
        an iterative method in its solution.  In that case,
        the residual vector Res = r - P z of the system
        should be made less than delta in weighted L2 norm,
        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
        Note: the error weight vector ewt can be obtained
        through a call to the routine ARKodeGetErrWeights.

 lr     is an input flag indicating whether MPrecSolve is to use
        the left preconditioner P1 or right preconditioner
        P2: lr = 1 means use P1, and lr = 2 means use P2.

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

 Returned value:
 The value to be returned by the MPrecSolve function is a flag
 indicating whether it was successful.  This value should be
   0 if successful,
   negative for an unrecoverable error (integration is halted).
---------------------------------------------------------------*/
typedef int (*ARKSpilsMassPrecSolveFn)(realtype t, N_Vector r, 
                                       N_Vector z, realtype delta, 
                                       int lr, void *user_data);
  
  
/*===============================================================
  ARKSPILS Exported functions
===============================================================*/

/*---------------------------------------------------------------
 Required inputs for the ARKSPILS linear solver interface:

 ARKSpilsSetLinearSolver specifies the iterative SUNLinearSolver 
 object that ARKode should use.  This is required if ARKode is 
 solving an implicit or IMEX IVP, and using an inexact Newton 
 solver.

 ARKSpilsSetMassLinearSolver specifies the iterative 
 SUNLinearSolver object that ARKode should use for mass-matrix 
 systems.  This is required if ARKode is solving a problem with 
 non-identity mass matrix and the user wishes to use an iterative 
 solver for these systems.

 NOTE: when solving an implicit or IMEX IVP with non-identity mass
 matrix and iterative linear solver, both the system and mass solvers
 must be iterative (i.e. you cannot combine a direct system 
 solver with an iterative mass matrix solver, etc.).

 The return value is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the ARKODE memory was NULL
    ARKSPILS_ILL_INPUT if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSpilsSetLinearSolver(void *arkode_mem, 
                                            SUNLinearSolver LS);

SUNDIALS_EXPORT int ARKSpilsSetMassLinearSolver(void *arkode_mem, 
                                                SUNLinearSolver LS,
                                                booleantype time_dep);

  
/*---------------------------------------------------------------
 Optional inputs to the ARKSPILS linear solver -- ALL of these 
 must be called AFTER the corresponding iterative linear solver 
 object (system matrix or mass matrix) has been attached to the 
 ARKode integrator).

 ARKSpilsSetEpsLin specifies the factor by which the tolerance on
                the nonlinear iteration is multiplied to get a
                tolerance on the linear iteration.
                Default value is 0.05.

 ARKSpilsSetMassEpsLin specifies the factor by which the tolerance
                on the nonlinear iteration is multiplied to get a
                tolerance on the mass matrix linear iteration.
                Default value is 0.05.

 ARKSpilsSetPreconditioner specifies the PrecSetup and PrecSolve 
                functions.  Default is NULL for both arguments 
                (no preconditioning)

 ARKSpilsSetMassPreconditioner specifies the mass matrix MPrecSetup
                and MPrecSolve functions.  Default is NULL for 
                both arguments (no preconditioning)

 ARKSpilsSetJacTimes specifies the jtsetup and jtimes functions. 
                Default is to use an internal finite difference 
		approximation routine (no setup).

 ARKSpilsSetMassTimes specifies the mtsetup and mtimes functions. 
                Note that there do not exist built-in finite-
                difference approximation routines for this, this 
                function MUST be called with non-NULL 'mtimes' 
                if ARKSpilsSetMassLinearSolver was called.

 The return value of ARKSpilsSet* is one of:
    ARKSPILS_SUCCESS      if successful
    ARKSPILS_MEM_NULL     if the arkode memory was NULL
    ARKSPILS_LMEM_NULL    if the linear solver memory was NULL
    ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was NULL
    ARKSPILS_ILL_INPUT    if an input has an illegal value
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSpilsSetEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int ARKSpilsSetMassEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int ARKSpilsSetPreconditioner(void *arkode_mem, 
                                              ARKSpilsPrecSetupFn psetup,
                                              ARKSpilsPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKSpilsSetMassPreconditioner(void *arkode_mem, 
                                                  ARKSpilsMassPrecSetupFn psetup,
                                                  ARKSpilsMassPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKSpilsSetJacTimes(void *arkode_mem, 
                                        ARKSpilsJacTimesSetupFn jtsetup,
                                        ARKSpilsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int ARKSpilsSetMassTimes(void *arkode_mem, 
                                         ARKSpilsMassTimesSetupFn msetup,
                                         ARKSpilsMassTimesVecFn mtimes,
                                         void *mtimes_data);

/*---------------------------------------------------------------
 Optional outputs from the ARKSPILS linear solver:

 ARKSpilsGetWorkSpace returns the real and integer workspace used
                by the SPILS module.

 ARKSpilsGetNumPrecEvals returns the number of preconditioner
                 evaluations, i.e. the number of calls made
                 to PrecSetup with jok==SUNFALSE.

 ARKSpilsGetNumPrecSolves returns the number of calls made to
                 PrecSolve.

 ARKSpilsGetNumLinIters returns the number of linear iterations.

 ARKSpilsGetNumConvFails returns the number of linear
                 convergence failures.

 ARKSpilsGetNumJTSetupEvals returns the number of calls to jtsetup.
 
 ARKSpilsGetNumJtimesEvals returns the number of calls to jtimes.

 ARKSpilsGetNumRhsEvals returns the number of calls to the user
                 f routine due to finite difference Jacobian
                 times vector evaluation.

 ARKSpilsGetLastFlag returns the last error flag set by any of
                 the ARKSPILS interface functions.

 ARKSpilsGetMassWorkSpace returns the real and integer workspace used
                by the mass matrix SPILS module.

 ARKSpilsGetNumMassPrecEvals returns the number of mass matrix 
                 preconditioner evaluations, i.e. the number of 
                 calls made to MPrecSetup.

 ARKSpilsGetNumMassPrecSolves returns the number of calls made to
                 MPrecSolve.

 ARKSpilsGetNumMassIters returns the number of mass matrix solver
                  iterations.

 ARKSpilsGetNumMassConvFails returns the number of mass matrix solver
                 convergence failures.

 ARKSpilsGetNumMtimesEvals returns the number of calls to mtimes.

 ARKSpilsGetLastMassFlag returns the last error flag set by any of
                 the ARKSPILS interface functions on the mass 
		 matrix solve.

 The return value of ARKSpilsGet* is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSpilsGetWorkSpace(void *arkode_mem, 
                                         long int *lenrwLS, 
                                         long int *leniwLS);
SUNDIALS_EXPORT int ARKSpilsGetNumPrecEvals(void *arkode_mem, 
                                            long int *npevals);
SUNDIALS_EXPORT int ARKSpilsGetNumPrecSolves(void *arkode_mem, 
                                             long int *npsolves);
SUNDIALS_EXPORT int ARKSpilsGetNumLinIters(void *arkode_mem, 
                                           long int *nliters);
SUNDIALS_EXPORT int ARKSpilsGetNumConvFails(void *arkode_mem, 
                                            long int *nlcfails);
SUNDIALS_EXPORT int ARKSpilsGetNumJTSetupEvals(void *arkode_mem,
                                               long int *njtsetups);
SUNDIALS_EXPORT int ARKSpilsGetNumJtimesEvals(void *arkode_mem, 
                                              long int *njvevals);
SUNDIALS_EXPORT int ARKSpilsGetNumRhsEvals(void *arkode_mem, 
                                           long int *nfevalsLS); 
SUNDIALS_EXPORT int ARKSpilsGetLastFlag(void *arkode_mem, 
                                        long int *flag);
  
SUNDIALS_EXPORT int ARKSpilsGetMassWorkSpace(void *arkode_mem, 
                                             long int *lenrwMLS, 
                                             long int *leniwMLS);
SUNDIALS_EXPORT int ARKSpilsGetNumMassPrecEvals(void *arkode_mem, 
                                                long int *nmpevals);
SUNDIALS_EXPORT int ARKSpilsGetNumMassPrecSolves(void *arkode_mem, 
                                                 long int *nmpsolves);
SUNDIALS_EXPORT int ARKSpilsGetNumMassIters(void *arkode_mem, 
                                            long int *nmiters);
SUNDIALS_EXPORT int ARKSpilsGetNumMassConvFails(void *arkode_mem, 
                                                long int *nmcfails);
SUNDIALS_EXPORT int ARKSpilsGetNumMtimesEvals(void *arkode_mem, 
                                              long int *nmvevals);
SUNDIALS_EXPORT int ARKSpilsGetLastMassFlag(void *arkode_mem, 
                                            long int *flag);


/*---------------------------------------------------------------
 The following function returns the name of the constant 
 associated with a ARKSPILS return flag
---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ARKSpilsGetReturnFlagName(long int flag);

#ifdef __cplusplus
}
#endif

#endif
