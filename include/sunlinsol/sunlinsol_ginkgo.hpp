/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * SUNLinearSolver interface to Ginkgo linear solvers
 * ---------------------------------------------------------------------------*/

#ifndef _SUNLINSOL_GINKGO_HPP
#define _SUNLINSOL_GINKGO_HPP

#include <cstring>
#include <memory>
#include <ginkgo/ginkgo.hpp>
#include <sundials/sundials_matrix.hpp>
#include <sundials/sundials_linearsolver.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

namespace sundials {
namespace ginkgo {

template<class GkoSolverType, class GkoMatrixType>
class LinearSolver;

inline SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S) { return SUNLINEARSOLVER_MATRIX_ITERATIVE; }

inline SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S) { return SUNLINEARSOLVER_GINKGO; }

inline int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S) { return SUNLS_SUCCESS; }

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
  auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  solver->setup(static_cast<Matrix<GkoMatrixType>*>(A->content));
  return SUNLS_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  solver->solve(b, x, tol);
  return SUNLS_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolFree_Ginkgo(SUNLinearSolver S)
{
  auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  delete solver; // NOLINT
  return SUNLS_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S)
{
  auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  return solver->numIters();
}

template<class GkoSolverType, class GkoMatrixType>
sunrealtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S)
{
  auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  return solver->resNorm();
}

// Custom gko::stop::Criterion that does the normal SUNDIALS stopping checks:
// 1. Was the absolute residual tolerance met?
// 2. Was the max iteration count reached?
class DefaultStop : public gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>
{
  friend class gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>;

public:
  GKO_CREATE_FACTORY_PARAMETERS(parameters, Factory)
  {
    sunrealtype GKO_FACTORY_PARAMETER_SCALAR(tolerance,
                                             SUN_UNIT_ROUNDOFF); // NOLINT(cppcoreguidelines-avoid-magic-numbers)
    gko::uint64 GKO_FACTORY_PARAMETER_SCALAR(max_iters, 50);
  };
  GKO_ENABLE_CRITERION_FACTORY(DefaultStop, parameters, Factory);
  GKO_ENABLE_BUILD_METHOD(Factory);

  gko::uint64 get_max_iters() const { return parameters_.max_iters; }

  sunrealtype get_tolerance() const { return parameters_.tolerance; }

protected:
  bool check_impl(gko::uint8 stoppingId, bool setFinalized, gko::array<gko::stopping_status>* stop_status,
                  bool* one_changed, const Updater&) override;

  explicit DefaultStop(std::shared_ptr<const gko::Executor> exec)
      : EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>(std::move(exec))
  {}

  explicit DefaultStop(const Factory* factory, const gko::stop::CriterionArgs& args)
      : EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>(factory->get_executor()), parameters_{
                                                                                                 factory->get_parameters()}
  {
    criteria_.push_back(gko::stop::AbsoluteResidualNorm<sunrealtype>::build() //
                            .with_tolerance(parameters_.tolerance)            //
                            .on(factory->get_executor())
                            ->generate(args));
    criteria_.push_back(gko::stop::Iteration::build()              //
                            .with_max_iters(parameters_.max_iters) //
                            .on(factory->get_executor())
                            ->generate(args));
  }

private:
  std::vector<std::unique_ptr<Criterion>> criteria_{};
};

inline bool DefaultStop::check_impl(gko::uint8 stoppingId, bool setFinalized,
                                    gko::array<gko::stopping_status>* stop_status, bool* one_changed,
                                    const Updater& updater)
{
  bool one_converged = false;
  gko::uint8 ids{1};
  *one_changed = false;
  for (auto& c : criteria_) {
    bool local_one_changed = false;
    one_converged |= c->check(ids, setFinalized, stop_status, &local_one_changed, updater);
    *one_changed |= local_one_changed;
    if (one_converged) {
      break;
    }
    ids++;
  }
  return one_converged;
}

template<class GkoSolverType, class GkoMatrixType>
class LinearSolver : public ConvertibleTo<SUNLinearSolver>
{
public:
  // Default constructor means the solver must be moved to
  LinearSolver() = default;

  LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, SUNContext sunctx)
      : gko_solver_factory_(gko_solver_factory), gko_solver_(nullptr),
        sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
        sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>()), res_norm_(sunrealtype{0.0})
  {
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    sunlinsol_->sunctx  = sunctx;

    sunlinsol_->ops->gettype    = SUNLinSolGetType_Ginkgo;
    sunlinsol_->ops->getid      = SUNLinSolGetID_Ginkgo;
    sunlinsol_->ops->initialize = SUNLinSolInitialize_Ginkgo;
    sunlinsol_->ops->setup      = SUNLinSolSetup_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->solve      = SUNLinSolSolve_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->numiters   = SUNLinSolNumIters_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->resnorm    = SUNLinSolResNorm_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->free       = SUNLinSolFree_Ginkgo<GkoSolverType, GkoMatrixType>;
  }

  // Don't allow copy constructor
  LinearSolver(const LinearSolver& that_solver) = delete;

  // Move constructor
  LinearSolver(LinearSolver&& that_solver) noexcept
      : gko_solver_factory_(std::move(that_solver.gko_solver_factory_)), gko_solver_(std::move(that_solver.gko_solver_)),
        sunlinsol_(std::move(that_solver.sunlinsol_)), sunlinsol_ops_(std::move(that_solver.sunlinsol_ops_)),
        iter_count_(that_solver.iter_count_), res_norm_(that_solver.res_norm_)
  {
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
  }

  // Don't allow copy assignment
  LinearSolver& operator=(const LinearSolver& rhs) = delete;

  // Move assignment
  LinearSolver& operator=(LinearSolver&& rhs)
  {
    gko_solver_factory_ = std::move(rhs.gko_solver_factory_);
    gko_solver_         = std::move(rhs.gko_solver_);
    sunlinsol_          = std::move(rhs.sunlinsol_);
    sunlinsol_ops_      = std::move(rhs.sunlinsol_ops_);
    iter_count_         = rhs.iter_count_;
    res_norm_           = rhs.res_norm_;
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    return *this;
  }

  ~LinearSolver() override = default;

  operator SUNLinearSolver() override { return sunlinsol_.get(); }
  operator SUNLinearSolver() const override { return sunlinsol_.get(); }
  SUNLinearSolver get() override { return sunlinsol_.get(); }
  SUNLinearSolver get() const override { return sunlinsol_.get(); }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gko_solver_factory_->get_executor(); }

  std::shared_ptr<typename GkoSolverType::Factory> gkofactory() { return gko_solver_factory_; }

  GkoSolverType* gkosolver() { return gko_solver_.get(); }

  int numIters() const { return iter_count_; }

  sunrealtype resNorm() const { return res_norm_; }

  GkoSolverType* setup(Matrix<GkoMatrixType>* A)
  {
    gko_solver_ = gko_solver_factory_->generate(A->gkomtx());
    return gko_solver_.get();
  }

  gko::LinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    auto logger{gko::share(gko::log::Convergence<sunrealtype>::create())};

    // Ginkgo provides a lot of options for stopping criterion,
    // so we make it possible to use it, but default to using
    // our normal iterative linear solver criterion.
    // If the criterion on the solver is of type DefaultStop,
    // then we will override the tolerance.
    auto crit{dynamic_cast<const DefaultStop::Factory*>(gkosolver()->get_stop_criterion_factory().get())};
    if (crit != nullptr) {
      auto new_crit = DefaultStop::build()     //
                          .with_tolerance(tol) //
                          .with_max_iters(crit->get_parameters().max_iters)
                          .on(gkoexec());
      new_crit->add_logger(logger);
      gkosolver()->set_stop_criterion_factory(std::move(new_crit));
    }

    gko::LinOp* result{nullptr};
    if (x != b) {
      auto x_vec = WrapVector(gkoexec(), x);
      auto b_vec = WrapVector(gkoexec(), b);

      // x = A^{-1} b
      result = gkosolver()->apply(b_vec.get(), x_vec.get());
    }
    else {
      auto x_vec = WrapVector(gkoexec(), x);

      // x = A^{-1} x
      result = gkosolver()->apply(x_vec.get(), x_vec.get());
    }

    iter_count_ += logger->get_num_iterations();
    gkoexec()->get_master()->copy_from(gko::lend(gkoexec()),
                                       1,
                                       gko::as<GkoDenseMat>(logger->get_residual_norm())->get_const_values(),
                                       &res_norm_);

    return result;
  }

private:
  std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory_;
  std::unique_ptr<GkoSolverType> gko_solver_;
  std::unique_ptr<_generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<_generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  int iter_count_{};
  sunrealtype res_norm_{};
};

} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGO_HPP
