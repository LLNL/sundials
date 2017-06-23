/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the header file for a generic linear solver package.
 * It defines the SUNLinearSolver structure (_generic_SUNLinearSolver)
 * which contains the following fields:
 *   - an implementation-dependent 'content' field which contains
 *     any internal data required by the solver
 *   - an 'ops' filed which contains a structure listing operations
 *     acting on/by such solvers
 *
 * Part I of this file contains enumeration constants for all 
 * SUNDIALS-defined linear solver types, as well as a generic type for 
 * user-supplied linear solver types.
 *
 * Part II of this file contains type declarations for the
 * _generic_SUNLinearSolver and _generic_SUNLinearSolver_Ops structures, 
 * as well as references to pointers to such structures 
 * (SUNLinearSolver).
 *
 * Part III of this file contains the prototypes for the linear solver
 * functions which operate on/by SUNLinearSolver objects.
 *
 * At a minimum, a particular implementation of a SUNLinearSolver must
 * do the following:
 *  - specify the 'content' field of SUNLinearSolver,
 *  - implement the operations on/by those SUNLinearSolver,
 *  - provide a constructor routine for new SUNLinearSolver objects
 *
 * Additionally, a SUNLinearSolver implementation may provide the 
 * following:
 *  - "Set" routines to control solver-specific parameters/options
 *  - "Get" routines to access solver-specific performance metrics
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINEARSOLVER_H
#define _SUNLINEARSOLVER_H
 
#include <sundials/sundials_types.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
  
/*
 * -----------------------------------------------------------------
 * I. Implemented SUNLinearSolver types: 
 *
 * These type names may be modified, but a at a minimum a client 
 * nonlinear solver and/or time integrator will want to know whether 
 * matrix/factorization information can be reused (hence the DIRECT 
 * and ITERATIVE types).
 * -----------------------------------------------------------------
 */
 
typedef enum {
  SUNLINEARSOLVER_DIRECT,
  SUNLINEARSOLVER_ITERATIVE,
  SUNLINEARSOLVER_CUSTOM
} SUNLinearSolver_Type;

  
/* 
 * -----------------------------------------------------------------
 * II. Generic definition of SUNLinearSolver 
 * -----------------------------------------------------------------
 */
 
/* Forward reference for pointer to SUNLinearSolver_Ops object */
typedef struct _generic_SUNLinearSolver_Ops *SUNLinearSolver_Ops;
 
/* Forward reference for pointer to SUNLinearSolver object */
typedef struct _generic_SUNLinearSolver *SUNLinearSolver;
 
/* Structure containing function pointers to linear solver operations */  
struct _generic_SUNLinearSolver_Ops {
  SUNLinearSolver_Type (*gettype)(SUNLinearSolver);
  int                  (*setatimes)(SUNLinearSolver, void*, 
                                    ATimesFn);
  int                  (*setpreconditioner)(SUNLinearSolver, void*, 
                                            PSetupFn, PSolveFn);
  int                  (*initialize)(SUNLinearSolver);
  int                  (*setup)(SUNLinearSolver, SUNMatrix, 
                                N_Vector, N_Vector, N_Vector);
  int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector, 
                                N_Vector, N_Vector, realtype);
  int                  (*performance)(SUNLinearSolver, int);
  long int             (*lastflag)(SUNLinearSolver);
  int                  (*free)(SUNLinearSolver);
};
 
/* A linear solver is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of linear solver
   operations corresponding to that implementation. */
struct _generic_SUNLinearSolver {
  void *content;
  struct _generic_SUNLinearSolver_Ops *ops;
};

  
/*
 * -----------------------------------------------------------------
 * III. Functions exported by SUNLinearSolver module
 * 
 * SUNLinSolGetType
 *   Returns an identifier for the linear solver type from 
 *   enumeration SUNLinearSolver_Type.
 *
 * SUNLinSolSetATimes (iterative methods only)
 *   Sets the function pointer for ATimes inside of an iterative 
 *   linear solver object.  This function should only be called by
 *   a main integrator, who will either provide this via 
 *   difference-quotients and vector operations, or by translating 
 *   between the generic PSetup and PSolve calls and the integrator-
 *   specific user-supplied routines.
 *
 * SUNLinSolSetPreconditioner (iterative methods only)
 *   Sets function pointers for PSetup and PSolve routines inside 
 *   of iterative linear solver objects.  This function should only 
 *   be called by a main integrator, who will provide translation 
 *   between the generic PSetup and PSolve calls and the integrator-
 *   specific user-supplied routines.
 *
 * SUNLinSolInitialize
 *   Performs linear solver initialization (assumes that all 
 *   solver-specific options have been set)
 *
 * SUNLinSolSetup
 *   Performs any linear solver setup needed, based on an updated
 *   system matrix A.  This may be called frequently (e.g. with a
 *   full Newton method) or infrequently (for a modified Newton 
 *   method), based on the type of integrator and/or nonlinear 
 *   solver requesting the solves.
 *
 * SUNLinSolSolve
 *   Solves a linear system A*x = b.  If the solver is scaled, it 
 *   uses the supplied scaling vector.  If the solver is iterative, 
 *   it attempts to solve to the specified tolerance (weighted RMS 
 *   norm).  If the solver is direct it ignores the input tolerance 
 *   and scaling vectors, and if the solver does not support scaling 
 *   then it should just use an RMS norm.
 *
 * SUNLinSolPerformance (optional; only required for IDA(s))
 *   Performs linear solver performance-related tasks: for perftask=0, 
 *   an initialization of performance variables is performed, while 
 *   for perftask=1, the performance is evaluated.
 *
 * SUNLinSolLastFlag
 *   Returns the last error flag encountered within the linear solver,
 *   allowing the user to investigate linear solver issues after 
 *   failed solves.
 *
 * SUNLinSolFree
 *   Frees memory allocated by the linear solver.
 *
 * ---------------------------------------------------------------
 *
 * The following table lists the linear solver functions used by 
 * different modules in SUNDIALS.  The symbols in the table have 
 * The following meaning:
 * M - called by the modified Newton nonlinear solver
 * I - called by the inexact Newton nonlinear solver
 * P - called by the Picard nonlinear solver
 * S - called by the integrator directly (and non-identity mass matrix)
 *
 *  LinearSolver
 *  Functions           CVODE(S)  ARKode     IDA(S)      KINSOL
 *  ------------------------------------------------------------
 *  GetType             M I P#    M I P#     M I P#      M I P
 *  SetATimes           I         I S+       I           I
 *  SetPreconditioner   I*        I* S*      I*          I*
 *  Initialize          M I P#    M I P# S   M I P#      M I P
 *  Setup               M I P#    M I P# S   M I P#      M I P      
 *  Solve               M I P#    M I P# S   M I P#      M I P
 *  Performance                              M I P#
 *  LastFlag^
 *  Free                M I P#    M I P# S   M I P#      M I P
 *  ------------------------------------------------------------
 * Notes: * -- only if user calls integrator-specific 
 *             preconditioner "set" routine
 *        + -- only called when using a non-identity mass matrix
 *             with an iterative linear solver
 *        # -- planned (currently only available for KINSOL)         
 *        ^ -- available for users to diagnose solver failures
 * ---------------------------------------------------------------
 */
  
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver S);

SUNDIALS_EXPORT int SUNLinSolSetATimes(SUNLinearSolver S, void* A_data,
                                       ATimesFn At);
  
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                               PSetupFn Pset, PSolveFn Psol);
  
SUNDIALS_EXPORT int SUNLinSolInitialize(SUNLinearSolver S);
  
SUNDIALS_EXPORT int SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A, N_Vector tmp1,
                                   N_Vector tmp2, N_Vector tmp3);
  
SUNDIALS_EXPORT int SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                   N_Vector b, N_Vector w, realtype tol);
  
SUNDIALS_EXPORT int SUNLinSolPerformance(SUNLinearSolver S, int perftask);
  
SUNDIALS_EXPORT long int SUNLinSolLastFlag(SUNLinearSolver S);
  
SUNDIALS_EXPORT int SUNLinSolFree(SUNLinearSolver S);
 
#ifdef __cplusplus
}
#endif
#endif
