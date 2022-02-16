#include <sunlinsol/sunlinsol_ginkgo.hpp>

using namespace sundials;

using GkoBatchDenseMat = gko::matrix::BatchDense<sunrealtype>;

#define GET_CONTENT(S)                                           \
  ((ginkgo::BlockLinearSolver<gko::solver::BatchCg<sunrealtype>, \
                              GkoBatchDenseMat>*)S->content)

SUNLinearSolver_Type SUNLinSolGetType_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

SUNLinearSolver_ID SUNLinSolGetID_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGO;
}

int SUNLinSolInitialize_GinkgoBlock(SUNLinearSolver S)
{
  return SUNMAT_SUCCESS;
}

int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S, N_Vector s1,
                                           N_Vector s2)
{
  auto solver = GET_CONTENT(S);
  solver->setScalingVectors(s1, s2);
  return SUNMAT_SUCCESS;
}

// Irrelavant since matrix-iterative.
// int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
//                                  booleantype onff);

int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = GET_CONTENT(S);
  solver->setup(static_cast<ginkgo::BlockMatrix<GkoBatchDenseMat>*>(A->content)
                    ->gkomtx());
  return SUNMAT_SUCCESS;
}

int SUNLinSolSolve_GinkgoBlock(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  auto solver = GET_CONTENT(S);

  // Ginkgo provides a lot of options for stopping criterion,
  // so we make it possible to use it, but default to using
  // our normal iterative linear solver criterion.
  if (!solver->useCustomStop())
  {
    solver->gkosolver()->set_stop_criterion_factory(
        gko::stop::AbsoluteResidualNorm<sunrealtype>::build().with_tolerance(
            tol));
  }

  solver->solve(b, x);
}

// We do not implement the following optional methods because it requires
// too much involvement in the solver building process. Users can get
// this info through ginkgo directly.
// int SUNLinSolNumIters_GinkgoBlock(SUNLinearSolver S);
// sunrealtype SUNLinSolResNorm_GinkgoBlock(SUNLinearSolver S);
// N_Vector SUNLinSolResid_GinkgoBlock(SUNLinearSolver S);
// sunindextype SUNLinSolLastFlag_GinkgoBlock(SUNLinearSolver S);
// int SUNLinSolSpace_GinkgoBlock(SUNLinearSolver S,
//                           long int *lenrwLS,
//                           long int *leniwLS);

int SUNLinSolFree_GinkgoBlock(SUNLinearSolver S) {}