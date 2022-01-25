
#include <memory>

#include <sundials/sundials_linearsolver.h>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

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
                                                     N_Vector s1,
                                                     N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
                                                booleantype onff);
SUNDIALS_EXPORT int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_Ginkgo(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_Ginkgo(SUNLinearSolver S,
                                         long int *lenrwLS,
                                         long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_Ginkgo(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

namespace sundials
{
namespace ginkgo
{

template<typename GkoSolverType, typename GkoMatType, typename GkoVecType>
class BaseLinearSolver
{

public:
  BaseLinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, bool use_custom_stop, SUNContext sunctx)
    : gkosolver_factory_(gko_solver_factory),
      sunlinsol_(SUNLinSolNewEmpty(sunctx)),
      use_custom_stop_(false)
  {
    sunlinsol_->content = this;

    sunlinsol_->ops->gettype           = SUNLinSolGetType_Ginkgo;
    sunlinsol_->ops->getid             = SUNLinSolGetID_Ginkgo;
    sunlinsol_->ops->setatimes         = SUNLinSolSetATimes_Ginkgo;
    sunlinsol_->ops->setpreconditioner = SUNLinSolSetPreconditioner_Ginkgo;
    sunlinsol_->ops->setscalingvectors = SUNLinSolSetScalingVectors_Ginkgo;
    sunlinsol_->ops->setzeroguess      = SUNLinSolSetZeroGuess_Ginkgo;
    sunlinsol_->ops->initialize        = SUNLinSolInitialize_Ginkgo;
    sunlinsol_->ops->setup             = SUNLinSolSetup_Ginkgo;
    sunlinsol_->ops->solve             = SUNLinSolSolve_Ginkgo;
    sunlinsol_->ops->numiters          = SUNLinSolNumIters_Ginkgo;
    sunlinsol_->ops->resnorm           = SUNLinSolResNorm_Ginkgo;
    sunlinsol_->ops->resid             = SUNLinSolResid_Ginkgo;
    sunlinsol_->ops->lastflag          = SUNLinSolLastFlag_Ginkgo;
    sunlinsol_->ops->space             = SUNLinSolSpace_Ginkgo;
    sunlinsol_->ops->free              = SUNLinSolFree_Ginkgo;
  }

  BaseLinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, SUNContext sunctx)
    : BaseLinearSolver(gko_solver_factory, gko_solver_factory, sunctx)
  {}

  std::shared_ptr<const gko::Executor> gkoexec() const { return gkosolver_->get_executor(); }
  std::shared_ptr<typename GkoSolverType::Factory> gkofactory() { return gkosolver_factory_; }
  GkoSolverType* gkosolver() { return gko::lend(gkosolver_); }
  bool useCustomStop() const { return use_custom_stop_; }

  GkoSolverType* setup(std::shared_ptr<GkoMatType> A)
  {
    gkosolver_ = gkosolver_factory_->generate(A);
    return gko::lend(gkosolver_);
  }

  gko::LinOp* solve(N_Vector b, N_Vector x)
  { throw std::runtime_error("sundials::ginkgo::LinearSolver::solve not implemented for specified gko::matrix type"); }

  void setScalingVectors(N_Vector left, N_Vector right)
  { throw std::runtime_error("sundials::ginkgo::LinearSolver::setScalingVectors not implemented for specified gko::matrix type"); }

protected:
  std::shared_ptr<typename GkoSolverType::Factory> gkosolver_factory_;
  std::unique_ptr<GkoSolverType> gkosolver_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  bool use_custom_stop_;

};

template<typename GkoSolverType, typename GkoMatType>
class LinearSolver : public BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>
{

public:
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::BaseLinearSolver;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkoexec;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkofactory;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkosolver;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::useCustomStop;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::setup;

  gko::LinOp* solve(N_Vector b, N_Vector x)
  {
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

};

template<typename GkoSolverType, typename GkoMatType>
class BatchLinearSolver : public BaseLinearSolver<GkoSolverType, GkoMatType, GkoBatchVecType>
{

public:
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::BaseLinearSolver;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkoexec;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkofactory;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::gkosolver;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::useCustomStop;
  using BaseLinearSolver<GkoSolverType, GkoMatType, GkoVecType>::setup;

  // gko::LinOp* solve(N_Vector b, N_Vector x)
  // {
  //   gko::LinOp* result;
  //   if (x != b)
  //   {
  //     auto x_vec = WrapVector(gkoexec(), x);
  //     auto b_vec = WrapVector(gkoexec(), b);

  //     // x = A^{-1} b
  //     result = gkosolver_->apply(gko::lend(b_vec), gko::lend(x_vec));
  //   }
  //   else
  //   {
  //     auto x_vec = WrapVector(gkoexec(), x);

  //     // x = A^{-1} x
  //     result = gkosolver_->apply(gko::lend(x_vec), gko::lend(x_vec));
  //   }

  //   return result;
  // }
};


}// namespace ginkgo
}// namespace sundials

#endif
