
#include <memory>
#include <sundials/sundials_linearsolver.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#ifndef _SUNLINSOL_GINKGO_HPP
#define _SUNLINSOL_GINKGO_HPP

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_Ginkgo(SUNLinearSolver S, void* A_data,
                                              SUNATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_Ginkgo(SUNLinearSolver S,
                                                      void* P_data,
                                                      SUNPSetupFn Pset,
                                                      SUNPSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S,
                                                      N_Vector s1, N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
                                                 booleantype onff);
SUNDIALS_EXPORT int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_Ginkgo(SUNLinearSolver S, long int* lenrwLS,
                                          long int* leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_Ginkgo(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

namespace sundials {
namespace ginkgo {

template<class GkoSolverType, class MatrixType> class LinearSolver;

SUNDIALS_EXPORT
SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

SUNDIALS_EXPORT
SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGO;
}

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S)
{
  return SUNMAT_SUCCESS;
}

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S,
                                                      N_Vector s1, N_Vector s2)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  return -1;
}

// Irrelavant since matrix-iterative.
// int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
//                                  booleantype onff);

template<typename LinearSolverType, typename MatrixType>
SUNDIALS_EXPORT int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  solver->setup(static_cast<MatrixType*>(A->content));
  return SUNLS_SUCCESS;
}

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  solver->solve(b, x, tol);
  return SUNLS_SUCCESS;
}

SUNDIALS_EXPORT
int SUNLinSolFree_Ginkgo(SUNLinearSolver S) { return SUNLS_SUCCESS; }

// Use SUNDIALS default stopping criteria.
std::unique_ptr<gko::stop::Combined::Factory> DefaultStop(
    std::shared_ptr<const gko::Executor> gko_exec,
    sunrealtype tol = SUN_UNIT_ROUNDOFF)
{
  auto tolerance_crit = gko::stop::AbsoluteResidualNorm<sunrealtype>::build() //
                            .with_tolerance(tol)                              //
                            .on(gko_exec);
  // TODO: In the unit tests, the max iters needed to be > 25 to get it to pass
  // auto iteration_crit = gko::stop::Iteration::build()         //
  //                           .with_max_iters(sunindextype{50}) //
  //                           .on(gko_exec);
  auto combined_crit = gko::stop::Combined::build()                  //
                           .with_criteria(std::move(tolerance_crit)) //
                           .on(gko_exec);

  // std::shared_ptr<const gko::log::Convergence<sunrealtype>> logger =
  // gko::log::Convergence<sunrealtype>::create(gko_exec);
  // combined_crit->add_logger(logger);

  return combined_crit;
}

template<class GkoSolverType, class MatrixType> class LinearSolver {
public:
  LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory,
               bool use_custom_criteria, SUNContext sunctx)
      : gko_solver_factory_(gko_solver_factory),
        sunlinsol_(SUNLinSolNewEmpty(sunctx)),
        use_custom_criteria_(use_custom_criteria)
  {
    // Attach function pointers for SUNLinearSolver
    using this_type = LinearSolver<GkoSolverType, MatrixType>;

    sunlinsol_->content      = this;
    sunlinsol_->ops->gettype = SUNLinSolGetType_Ginkgo;
    sunlinsol_->ops->getid   = SUNLinSolGetID_Ginkgo;
    sunlinsol_->ops->setscalingvectors =
        SUNLinSolSetScalingVectors_Ginkgo<this_type>;
    sunlinsol_->ops->initialize = SUNLinSolInitialize_Ginkgo<this_type>;
    sunlinsol_->ops->setup      = SUNLinSolSetup_Ginkgo<this_type, MatrixType>;
    sunlinsol_->ops->solve      = SUNLinSolSolve_Ginkgo<this_type>;
    // sunlinsol_->ops->lastflag          = SUNLinSolLastFlag_Ginkgo;
    // sunlinsol_->ops->space             = SUNLinSolSpace_Ginkgo;
    sunlinsol_->ops->free = SUNLinSolFree_Ginkgo;
  }

  LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory,
               SUNContext sunctx)
      : LinearSolver(gko_solver_factory, false, sunctx)
  {}

  std::shared_ptr<const gko::Executor> gkoexec() const
  {
    return gko_solver_->get_executor();
  }

  std::shared_ptr<typename GkoSolverType::Factory> gkofactory()
  {
    return gko_solver_factory_;
  }

  GkoSolverType* gkosolver() { return gko::lend(gko_solver_); }

  bool useCustomCriteria() const { return use_custom_criteria_; }

  operator SUNLinearSolver() { return sunlinsol_.get(); }

  operator SUNLinearSolver() const { return sunlinsol_.get(); }

  GkoSolverType* setup(MatrixType* A)
  {
    gko_solver_ = gko_solver_factory_->generate(A->gkomtx());
    return gko::lend(gko_solver_);
  }

  gko::LinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    // Ginkgo provides a lot of options for stopping criterion,
    // so we make it possible to use it, but default to using
    // our normal iterative linear solver criterion.
    if (!useCustomCriteria())
    {
      gkosolver()->set_stop_criterion_factory(
          gko::share(DefaultStop(gkoexec(), tol)));
    }

    gko::LinOp* result;
    if (x != b)
    {
      auto x_vec = WrapVector(gkoexec(), x);
      auto b_vec = WrapVector(gkoexec(), b);

      // x = A^{-1} b
      result = gkosolver()->apply(gko::lend(b_vec), gko::lend(x_vec));
    }
    else
    {
      auto x_vec = WrapVector(gkoexec(), x);

      // x = A^{-1} x
      result = gkosolver()->apply(gko::lend(x_vec), gko::lend(x_vec));
    }

    return result;
  }

protected:
  std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory_;
  std::unique_ptr<GkoSolverType> gko_solver_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  bool use_custom_criteria_;
};

// template<typename GkoSolverType, typename GkoMatType>
// class BlockLinearSolver : public BaseLinearSolver<GkoSolverType, GkoMatType,
// GkoBatchVecType>
// {

// public:
//   using BaseLinearSolver<GkoSolverType, GkoMatType,
//   GkoVecType>::BaseLinearSolver; using BaseLinearSolver<GkoSolverType,
//   GkoMatType, GkoVecType>::gkoexec; using BaseLinearSolver<GkoSolverType,
//   GkoMatType, GkoVecType>::gkofactory; using BaseLinearSolver<GkoSolverType,
//   GkoMatType, GkoVecType>::gkosolver; using BaseLinearSolver<GkoSolverType,
//   GkoMatType, GkoVecType>::useCustomCriteria; using
//   BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::setup;

//   // gko::LinOp* solve(N_Vector b, N_Vector x)
//   // {
//   //   gko::LinOp* result;
//   //   if (x != b)
//   //   {
//   //     auto x_vec = WrapVector(gkoexec(), x);
//   //     auto b_vec = WrapVector(gkoexec(), b);

//   //     // x = A^{-1} b
//   //     result = gko_solver_->apply(gko::lend(b_vec), gko::lend(x_vec));
//   //   }
//   //   else
//   //   {
//   //     auto x_vec = WrapVector(gkoexec(), x);

//   //     // x = A^{-1} x
//   //     result = gko_solver_->apply(gko::lend(x_vec), gko::lend(x_vec));
//   //   }

//   //   return result;
//   // }
// };

} // namespace ginkgo
} // namespace sundials

#endif
