/* -----------------------------------------------------------------------------
 * Programmer(s): David Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for a generic SUNNONLINEARSOLVER package. It
 * contains the implementation of the SUNNonlinearSolver operations listed in
 * sundials_nonlinearsolver.h
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/sundials_nonlinearsolver.h>

/* -----------------------------------------------------------------------------
 * Functions in the 'ops' structure
 * ---------------------------------------------------------------------------*/

/*
 * core functions
 */

SUNNonlinearSolver_Type SUNNonlinSolGetType(SUNNonlinearSolver NLS)
{
  return(NLS->ops->gettype(NLS));
}

int SUNNonlinSolInit(SUNNonlinearSolver NLS, N_Vector tmpl)
{
  return((int) NLS->ops->init(NLS, tmpl));
}

int SUNNonlinSolSetup(SUNNonlinearSolver NLS, N_Vector y, void* mem)
{
  return((int) NLS->ops->setup(NLS, y, mem));
}

int SUNNonlinSolSolve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                      N_Vector w, realtype tol, void* mem)
{
  return((int) NLS->ops->solve(NLS, y0, y, w, tol, mem));
}

int SUNNonlinSolFree(SUNNonlinearSolver NLS)
{
  if (NLS == NULL) return(0);
  return(NLS->ops->free(NLS));
}

/*
 * set functions
 */

int SUNNonlinSolSetSysFn(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  return((int) NLS->ops->setsysfn(NLS, SysFn));
}

int SUNNonlinSolSetLSetupFn(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn LSetupFn)
{
  return((int) NLS->ops->setlsetupfn(NLS, LSetupFn));
}

int SUNNonlinSolSetLSolveFn(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn LSolveFn)
{
  return((int) NLS->ops->setlsolvefn(NLS, LSolveFn));
}

int SUNNonlinSolSetConvTestFn(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn)
{
  return((int) NLS->ops->setctestfn(NLS, CTestFn));
}

int SUNNonlinSolSetMaxIters(SUNNonlinearSolver NLS, int maxiters)
{
  return((int) NLS->ops->setmaxiters(NLS, maxiters));
}

/*
 * get functions
 */

int SUNNonlinSolGetNumIters(SUNNonlinearSolver NLS, long int *niters)
{
  return((int) NLS->ops->getnumiters(NLS, niters));
}
