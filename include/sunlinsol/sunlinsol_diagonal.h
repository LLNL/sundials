/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the header file for the diagonal implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the diagonal 
 * implementation of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNDiagonalLinearSolver as well as implementation-specific 
 * prototypes for various useful matrix operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_DIAGONAL_H
#define _SUNLINSOL_DIAGONAL_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_diagonal.h>


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Diagonal implementation of SUNLinearSolver
 *
 * The diagonal implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     last_flag -- last error return flag from internal setup/solve
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_Diagonal {
  long int last_flag;
};

typedef struct _SUNLinearSolverContent_Diagonal *SUNLinearSolverContent_Diagonal;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_diagonal
 * 
 * CONSTRUCTORS:
 *    SUNDiagonalLinearSolver creates and allocates memory for a 
 * 	 diagonal matrix solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNDiagonalLinearSolver(N_Vector y,
                                                        SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * diagonal implementations of various linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_Diagonal(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_Diagonal(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_Diagonal(SUNLinearSolver S, void* A_data,
                                                ATSetupFn ATSetup, ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_Diagonal(SUNLinearSolver S,
                                                        void* P_data,
                                                        PSetupFn Pset,
                                                        PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_Diagonal(SUNLinearSolver S,
                                                        N_Vector s1,
                                                        N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetup_Diagonal(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Diagonal(SUNLinearSolver S, SUNMatrix A,
                                            N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_Diagonal(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_Diagonal(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_Diagonal(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_Diagonal(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif

