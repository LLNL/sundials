/*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
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
 * This is the header file contains the prototypes for functions to 
 * test SUNLinearSolver module implementations. 
 * -----------------------------------------------------------------
 */

#include <math.h>

/* define constatnts */
#define ZERO     RCONST(0.0)
#define ONE      RCONST(1.0)

/* NAN and floating point "equality" check, failure update macro */
#if __STDC_VERSION__ >= 199901L
#define FNEQ(a,b,tol) (isnan(a) ? 1 : ( SUNRabs((a)-(b)) > tol ))
#else
#define FNEQ(a,b,tol) (( SUNRabs((a)-(b)) > tol ))
#endif


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  /* Forward declarations for implementation specific utility functions */
  int check_vector(N_Vector X, N_Vector Y, realtype tol);
   
  /* Test function declarations */
  int Test_SUNLinSolGetType(SUNLinearSolver S, SUNLinearSolver_Type sunid, int myid);
  int Test_SUNLinSolLastFlag(SUNLinearSolver S, int myid);
  int Test_SUNLinSolNumIters(SUNLinearSolver S, int myid);
  int Test_SUNLinSolSetATimes(SUNLinearSolver S, booleantype shouldpass, int myid);
  int Test_SUNLinSolSetPreconditioner(SUNLinearSolver S, booleantype shouldpass, int myid);
  int Test_SUNLinSolInitialize(SUNLinearSolver S, int myid);
  int Test_SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A, int myid);
  int Test_SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                          N_Vector b, N_Vector w, realtype tol, int myid);

  /* Timing function */
  void SetTiming(int onoff);

  /* Dummy integrator-specific function pointers */
  int ATSetup(void *A_data);
  int ATimes(void *A_data, N_Vector v, N_Vector z);
  int PSetup(void *P_data);
  int PSolve(void *P_data, N_Vector r, N_Vector z, realtype tol, int lr);
  
#ifdef __cplusplus
}
#endif
