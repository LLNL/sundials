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
 * This is the implementation file for a generic SUNLINEARSOLVER 
 * package.  It contains the implementation of the SUNLinearSolver
 * operations listed in sundials_linearsolver.h
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include <sundials/sundials_linearsolver.h>

/*
 * -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver S)
{
  SUNLinearSolver_Type type;
  type = S->ops->gettype(S);
  return(type);
}

int SUNLinSolSetATimes(SUNLinearSolver S, void* A_data,
                       ATSetupFn ATSetup, ATimesFn ATimes)
{
  return ((int) S->ops->setatimes(S, A_data, ATSetup, ATimes));
}

  
int SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                               PSetupFn Pset, PSolveFn Psol)
{
  return ((int) S->ops->setpreconditioner(S, P_data, Pset, Psol));
}
  
int SUNLinSolSetScalingVectors(SUNLinearSolver S,
                               N_Vector s1, N_Vector s2)
{
  return ((int) S->ops->setscalingvectors(S, s1, s2));
}
  
int SUNLinSolInitialize(SUNLinearSolver S)
{
  return ((int) S->ops->initialize(S));
}
  
int SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A)
{
  return ((int) S->ops->setup(S, A));
}

int SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                   N_Vector b, realtype tol)
{
  return ((int) S->ops->solve(S, A, x, b, tol));
}
  
int SUNLinSolNumIters(SUNLinearSolver S)
{
  return ((int) S->ops->numiters(S));
}

realtype SUNLinSolResNorm(SUNLinearSolver S)
{
  return ((realtype) S->ops->resnorm(S));
}

int SUNLinSolNumPSolves(SUNLinearSolver S)
{
  return ((int) S->ops->numpsolves(S));
}

long int SUNLinSolLastFlag(SUNLinearSolver S)
{
  return ((long int) S->ops->lastflag(S));
}

int SUNLinSolFree(SUNLinearSolver S)
{
  if (S==NULL) return 0;
  S->ops->free(S);
  return 0;
}

