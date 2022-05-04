
#include <cstring>
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
SUNDIALS_EXPORT int SUNLinSolSetATimes_Ginkgo(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_Ginkgo(SUNLinearSolver S, void* P_data, SUNPSetupFn Pset,
                                                      SUNPSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S, N_Vector s1, N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S, booleantype onff);
SUNDIALS_EXPORT int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_Ginkgo(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_Ginkgo(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

namespace sundials {
namespace ginkgo {

template<class GkoSolverType, class MatrixType>
class LinearSolver;

SUNDIALS_EXPORT
SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S) { return SUNLINEARSOLVER_MATRIX_ITERATIVE; }

SUNDIALS_EXPORT
SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S) { return SUNLINEARSOLVER_GINKGO; }

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S)
{
  return SUNMAT_SUCCESS;
}

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  // auto solver = static_cast<LinearSolverType*>(S->content);

  return -1;
}

template<typename LinearSolverType, typename MatrixType>
SUNDIALS_EXPORT int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  solver->setup(static_cast<MatrixType*>(A->content));
  return SUNLS_SUCCESS;
}

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  solver->solve(b, x, tol);
  return SUNLS_SUCCESS;
}

SUNDIALS_EXPORT
int SUNLinSolFree_Ginkgo(SUNLinearSolver S) { return SUNLS_SUCCESS; }

template<typename LinearSolverType>
SUNDIALS_EXPORT int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  return solver->numIters();
}

template<typename LinearSolverType>
SUNDIALS_EXPORT sunrealtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S)
{
  auto solver = static_cast<LinearSolverType*>(S->content);
  return solver->resNorm();
}

// Custom gko::stop::Criterion that does the normal SUNDIALS stopping checks:
// 1. Was the absolute residual tolerance met?
// 2. Was the max iteration count reached?
class DefaultStop : public gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion> {
  friend class gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>;

public:
  GKO_CREATE_FACTORY_PARAMETERS(parameters, Factory)
  {
    sunrealtype GKO_FACTORY_PARAMETER_SCALAR(tolerance, SUN_UNIT_ROUNDOFF);
    sunrealtype GKO_FACTORY_PARAMETER_SCALAR(max_iters, 50);
  };
  GKO_ENABLE_CRITERION_FACTORY(DefaultStop, parameters, Factory);
  GKO_ENABLE_BUILD_METHOD(Factory);

  int get_max_iters() const { return parameters_.max_iters; }

  sunrealtype get_tolerance() const { return parameters_.tolerance; }

protected:
  bool check_impl(gko::uint8 stoppingId, bool setFinalized, gko::Array<gko::stopping_status>* stop_status,
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

bool DefaultStop::check_impl(gko::uint8 stoppingId, bool setFinalized, gko::Array<gko::stopping_status>* stop_status,
                             bool* one_changed, const Updater& updater)
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

class LinearSolverInterface {
public:
  virtual gko::LinOp* solve(N_Vector b, N_Vector x, sunrealtype tol) = 0;
  virtual std::shared_ptr<const gko::Executor> gkoexec() const       = 0;
  virtual int numIters() const                                       = 0;
  virtual sunrealtype resNorm() const                                = 0;
  virtual bool logResNorm() const                                    = 0;
  virtual bool logResNorm(bool onoff)                                = 0;
  virtual SUNLinearSolver get()                                      = 0;
  virtual SUNLinearSolver get() const                                = 0;
  virtual operator SUNLinearSolver()                                 = 0;
  virtual operator SUNLinearSolver() const                           = 0;
};

template<class GkoSolverType, class MatrixType>
class LinearSolver : public LinearSolverInterface {
public:
  LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, SUNContext sunctx)
      : gko_solver_factory_(gko_solver_factory), sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
        sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>()), iter_count_(0), res_norm_(sunrealtype{0.0}),
        log_res_norm_(false)
  {
    // Attach function pointers for SUNLinearSolver
    using this_type = LinearSolver<GkoSolverType, MatrixType>;

    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    sunlinsol_->sunctx  = sunctx;

    std::memset(sunlinsol_->ops, 0, sizeof(_generic_SUNLinearSolver_Ops));
    sunlinsol_->ops->gettype           = SUNLinSolGetType_Ginkgo;
    sunlinsol_->ops->getid             = SUNLinSolGetID_Ginkgo;
    sunlinsol_->ops->setscalingvectors = SUNLinSolSetScalingVectors_Ginkgo<this_type>;
    sunlinsol_->ops->initialize        = SUNLinSolInitialize_Ginkgo<this_type>;
    sunlinsol_->ops->setup             = SUNLinSolSetup_Ginkgo<this_type, MatrixType>;
    sunlinsol_->ops->solve             = SUNLinSolSolve_Ginkgo<this_type>;
    sunlinsol_->ops->numiters          = SUNLinSolNumIters_Ginkgo<this_type>;
    sunlinsol_->ops->resnorm           = SUNLinSolResNorm_Ginkgo<this_type>;
    sunlinsol_->ops->free              = SUNLinSolFree_Ginkgo;
  }

  std::shared_ptr<const gko::Executor> gkoexec() const { return gko_solver_->get_executor(); }

  std::shared_ptr<typename GkoSolverType::Factory> gkofactory() { return gko_solver_factory_; }

  GkoSolverType* gkosolver() { return gko::lend(gko_solver_); }

  int numIters() const { return iter_count_; }

  sunrealtype resNorm() const { return res_norm_; }

  bool logResNorm() const { return log_res_norm_; }

  bool logResNorm(bool onoff)
  {
    log_res_norm_ = onoff;
    return log_res_norm_;
  }

  operator SUNLinearSolver() { return sunlinsol_.get(); }

  operator SUNLinearSolver() const { return sunlinsol_.get(); }

  SUNLinearSolver get() { return sunlinsol_.get(); }

  SUNLinearSolver get() const { return sunlinsol_.get(); }

  GkoSolverType* setup(MatrixType* A)
  {
    gko_solver_ = gko_solver_factory_->generate(A->gkomtx());
    return gko::lend(gko_solver_);
  }

  gko::LinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    auto logger = gko::share(gko::log::Convergence<sunrealtype>::create(gkoexec()));

    // Ginkgo provides a lot of options for stopping criterion,
    // so we make it possible to use it, but default to using
    // our normal iterative linear solver criterion.
    // If the criterion on the solver is of type DefaultStop,
    // then we will override the tolerance.
    auto crit = dynamic_cast<const DefaultStop::Factory*>(gkosolver()->get_stop_criterion_factory().get());
    if (crit != nullptr) {
      auto new_crit = DefaultStop::build()     //
                          .with_tolerance(tol) //
                          .with_max_iters(crit->get_parameters().max_iters)
                          .on(gkoexec());
      new_crit->add_logger(logger);
      gkosolver()->set_stop_criterion_factory(gko::share(new_crit));
    }

    gko::LinOp* result;
    if (x != b) {
      auto x_vec = WrapVector(gkoexec(), x);
      auto b_vec = WrapVector(gkoexec(), b);

      // x = A^{-1} b
      result = gkosolver()->apply(b_vec.get(), x_vec.get());
    } else {
      auto x_vec = WrapVector(gkoexec(), x);

      // x = A^{-1} x
      result = gkosolver()->apply(x_vec.get(), x_vec.get());
    }

    iter_count_ += logger->get_num_iterations();
    if (logResNorm())
      res_norm_ = gko::as<GkoDenseMat>(logger->get_residual_norm())->at(0);

    return result;
  }

private:
  std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory_;
  std::unique_ptr<GkoSolverType> gko_solver_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<struct _generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  int iter_count_;
  sunrealtype res_norm_;
  bool log_res_norm_;
};

} // namespace ginkgo
} // namespace sundials

#endif
