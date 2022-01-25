#include <sunlinsol/sunlinsol_ginkgo.hpp>

using namespace sundials;

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;

#define GET_CONTENT(S) ((ginkgo::LinearSolver<gko::solver::Cg<sunrealtype>, GkoDenseMat> *) S->content)

SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGO;
}

int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S)
{
  return SUNMAT_SUCCESS;
}

// We dont need this operation since it is matrix-iterative.
// int SUNLinSolSetATimes_Ginkgo(SUNLinearSolver S, void* A_data,
//                               SUNATimesFn ATimes);

// We dont provide this operation since the Ginkgo
// preconditioner interface should be used instead.
// int SUNLinSolSetPreconditioner_Ginkgo(SUNLinearSolver S,
//                                       void* P_data,
//                                       SUNPSetupFn Pset,
//                                       SUNPSolveFn Psol);

// int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S,
//                                            N_Vector s1,
//                                            N_Vector s2);

int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S,
                                           N_Vector s1,
                                           N_Vector s2)
{
  auto solver = GET_CONTENT(S);

}

// Irrelavant since matrix-iterative.
// int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
//                                  booleantype onff);

int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = GET_CONTENT(S);
  solver->setup(static_cast<ginkgo::Matrix<GkoDenseMat>*>(A->content)->gkomtx());
  return SUNMAT_SUCCESS;
}

int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A,
                          N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver = GET_CONTENT(S);

  // Ginkgo provides a lot of options for stopping criterion,
  // so we make it possible to use it, but default to using
  // our normal iterative linear solver criterion.
  if (!solver->useCustomStop())
  {
    auto new_crit = gko::stop::AbsoluteResidualNorm<sunrealtype>::build()
                        .with_tolerance(tol)
                        .on(solver->gkoexec());
    solver->gkosolver()->set_stop_criterion_factory(gko::share(new_crit));
  }

  solver->solve(b, x);
}

// We do not implement the following optional methods because it requires
// too much involvement in the solver building process. Users can get
// this info through ginkgo directly.
// int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S);
// sunrealtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S);
// N_Vector SUNLinSolResid_Ginkgo(SUNLinearSolver S);
// sunindextype SUNLinSolLastFlag_Ginkgo(SUNLinearSolver S);
// int SUNLinSolSpace_Ginkgo(SUNLinearSolver S,
//                           long int *lenrwLS,
//                           long int *leniwLS);

int SUNLinSolFree_Ginkgo(SUNLinearSolver S)
{

}