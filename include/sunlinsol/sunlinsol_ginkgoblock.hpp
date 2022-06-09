
#include <cstring>
#include <memory>
#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>
#include <sundials/sundials_logger.h>

#ifndef _SUNLINSOL_GINKGOBLOCK_HPP
#define _SUNLINSOL_GINKGOBLOCK_HPP

namespace sundials {
namespace ginkgo {
namespace {

template<class GkoBatchSolverType, class BlockMatrixType>
class BlockLinearSolver;

SUNDIALS_EXPORT
SUNLinearSolver_Type SUNLinSolGetType_GinkgoBlock(SUNLinearSolver S) { return SUNLINEARSOLVER_MATRIX_ITERATIVE; }

SUNDIALS_EXPORT
SUNLinearSolver_ID SUNLinSolGetID_GinkgoBlock(SUNLinearSolver S) { return SUNLINEARSOLVER_GINKGOBLOCK; }

template<typename GkoBatchLinearSolverType>
SUNDIALS_EXPORT int SUNLinSolInitialize_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLS_SUCCESS;
}

template<typename GkoBatchLinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  if (s1 == NULL || s2 == NULL) {
    return SUNLS_ILL_INPUT;
  }
  solver->setScalingVectors(s1, s2);
  return SUNLS_SUCCESS;
}

template<typename GkoBatchLinearSolverType, typename BlockMatrixType>
SUNDIALS_EXPORT int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->setup(static_cast<BlockMatrixType*>(A->content));
}

template<typename GkoBatchLinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSolve_GinkgoBlock(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  if (solver->solve(b, x, tol)) { return SUNLS_SUCCESS; }
  else                          { return SUNLS_RES_REDUCED; } // TODO: how can we know if it was reduced?
}

SUNDIALS_EXPORT
int SUNLinSolFree_GinkgoBlock(SUNLinearSolver S) { return SUNLS_SUCCESS; }

template<typename GkoBatchLinearSolverType>
SUNDIALS_EXPORT int SUNLinSolNumIters_GinkgoBlock(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return std::round(solver->sumAvgNumIters());
}

inline std::unique_ptr<gko::matrix::BatchDiagonal<sunrealtype>> WrapBatchDiagMatrix(
    std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks, N_Vector x)
{
  sunrealtype* x_arr          = (x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x);
  const sunindextype xvec_len = N_VGetLength(x);
  auto batch_xvec_size        = gko::batch_dim<>(num_blocks, gko::dim<2>(xvec_len / num_blocks, xvec_len / num_blocks));
  auto xvec_view              = gko::Array<sunrealtype>::view(gko_exec, xvec_len, x_arr);
  return gko::matrix::BatchDiagonal<sunrealtype>::create(gko_exec, batch_xvec_size, std::move(xvec_view));
}

template<class GkoBatchSolverType, class BlockMatrixType>
class BlockLinearSolver : public SUNLinearSolverView {
public:
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    std::shared_ptr<gko::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_blocks,
                    SUNContext sunctx)
      : gko_exec_(gko_exec), tolerance_type_(tolerance_type), precon_factory_(precon_factory), max_iters_(max_iters),
        sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
        sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>()), num_blocks_(num_blocks),
        left_scale_vec_(nullptr), right_scale_vec_(nullptr), matrix_(nullptr),
        log_res_norm_(false), avg_iter_count_(0.0), sum_of_avg_iters_(0.0), stddev_iter_count_(0.0)
  {
    initSUNLinSol(sunctx);
  }

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, nullptr,
                          default_max_iters_, num_blocks, sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, nullptr, default_max_iters_,
                          num_blocks, sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::BatchLinOpFactory> precon_factory,
                    sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, precon_factory, default_max_iters_,
                          num_blocks, sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, int max_iters, sunindextype num_blocks,
                    SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, nullptr,
                          max_iters, num_blocks, sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    int max_iters, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, nullptr, max_iters, num_blocks,
                          sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::BatchLinOpFactory> precon_factory,
                    int max_iters, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, precon_factory, max_iters, num_blocks,
                          sunctx)
  {}

  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    std::shared_ptr<gko::BatchLinOpFactory> precon_factory, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, precon_factory, default_max_iters_, num_blocks, sunctx)
  {}

  std::shared_ptr<const gko::Executor> gkoExec() const { return gko_exec_; }

  SUNLinearSolver get() override { return sunlinsol_.get(); }

  SUNLinearSolver get() const override { return sunlinsol_.get(); }

  operator SUNLinearSolver() override { return sunlinsol_.get(); }

  operator SUNLinearSolver() const override { return sunlinsol_.get(); }

  SUNLogger sunLogger() {
    SUNContext_GetLogger(sunlinsol_->sunctx, &sunlogger_);
    return sunlogger_;
  }

  bool logResNorm() const { return log_res_norm_; }

  bool logResNorm(bool onoff)
  {
    log_res_norm_ = onoff;
    return log_res_norm_;
  }

  sunrealtype avgNumIters() const { return avg_iter_count_; }

  sunrealtype stddevNumIters() const { return stddev_iter_count_; }

  sunrealtype sumAvgNumIters () const { return sum_of_avg_iters_; }

  void setScalingVectors(N_Vector s1, N_Vector s2)
  {
    left_scale_vec_  = gko::share(WrapBatchDiagMatrix(gkoExec(), num_blocks_, s1));
    right_scale_vec_ = gko::share(WrapBatchDiagMatrix(gkoExec(), num_blocks_, s2));
  }

  int setup(BlockMatrixType* A)
  {
    if (num_blocks_ != A->numBlocks()) {
      return SUNLS_ILL_INPUT;
    }
    matrix_ = A;
    return SUNLS_SUCCESS;
  }

  gko::BatchLinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    auto solver_factory = GkoBatchSolverType::build()                     //
                              .with_max_iterations(max_iters_)            //
                              .with_residual_tol(tol)                     //
                              .with_tolerance_type(tolerance_type_)       //
                              .with_preconditioner(precon_factory_)       //
                              .with_left_scaling_op(left_scale_vec_)      //
                              .with_right_scaling_op(right_scale_vec_)    //
                              .on(gkoExec());
    auto solver = solver_factory->generate(matrix_->gkomtx());
    // if (!logger_) {
    //   logger_  = gko::share(gko::log::BatchConvergence<sunrealtype>::create(gkoExec()));
    //   solver->add_logger(logger_);
    // }
    auto logger_ = gko::share(gko::log::BatchConvergence<sunrealtype>::create(gkoExec()));
    solver->add_logger(logger_);

    gko::BatchLinOp* result = nullptr;
    if (x != b) {
      auto x_vec = WrapBatchVector(gkoExec(), num_blocks_, x);
      auto b_vec = WrapBatchVector(gkoExec(), num_blocks_, b);

      // x = A'^{-1} diag(left) b
      result = solver->apply(b_vec.get(), x_vec.get());
    } else {
      auto x_vec = WrapBatchVector(gkoExec(), num_blocks_, x);

      // x = A^'{-1} diag(right) x
      result = solver->apply(x_vec.get(), x_vec.get());
    }

    solver->remove_logger(logger_.get());

    auto num_iters = logger_->get_num_iterations();
    num_iters.set_executor(gkoExec()->get_master()); // move data to host

    // Check if any batch entry reached the max number of iterations.
    // While we check, compute the average number of iterations across all batch entries.
    bool at_least_one_maxed_out = false;
    auto iter_count = num_iters.get_const_data();
    auto max_iter_count = 0;
    auto min_iter_count = max_iters_;
    avg_iter_count_ = 0.0;
    for (int i = 0; i < num_blocks_; i++) {
      avg_iter_count_ += iter_count[i];
      max_iter_count = std::max(max_iter_count, iter_count[i]);
      min_iter_count = std::min(min_iter_count, iter_count[i]);
      if (iter_count[i] >= max_iters_) {
        at_least_one_maxed_out = true;
      }
    }
    avg_iter_count_ /= num_blocks_;
    sum_of_avg_iters_ += avg_iter_count_;

    // Compute the std. dev. in iteration count across all batch enties.
    stddev_iter_count_ = 0.0;
    for (int i = 0; i < num_blocks_; i++) {
      stddev_iter_count_ += std::pow(std::abs(iter_count[i] - avg_iter_count_), 2);
    }
    stddev_iter_count_ = std::sqrt(stddev_iter_count_ / num_blocks_);

    if (logResNorm()) {
      auto res_norm = logger_->get_residual_norm();
      std::cout << ">>>>>>> ginkgo res_norms:\n";
      for (auto& mat : res_norm->unbatch()) {
        std::cout << mat->at(0) << ", ";
      }
      std::cout << "\n";
    }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(sunLogger(), SUN_LOGLEVEL_INFO,
      "sundials::ginkgo::BlockLinearSolver::solve", "end-solve",
      "avg. iter count = %.16g, stddev. iter count = %.16g, max iter count = %d, min iter count = %d", avg_iter_count_, stddev_iter_count_, max_iter_count, min_iter_count);
#endif

    if (at_least_one_maxed_out) {
      return nullptr; // return no result to indicate we did not converge
    } else {
      return result;
    }
  }

private:
  std::shared_ptr<const gko::Executor> gko_exec_;
  gko::stop::batch::ToleranceType tolerance_type_;
  std::shared_ptr<gko::BatchLinOpFactory> precon_factory_;
  int max_iters_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<struct _generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  sunindextype num_blocks_;
  std::shared_ptr<gko::matrix::BatchDiagonal<sunrealtype>> left_scale_vec_;
  std::shared_ptr<gko::matrix::BatchDiagonal<sunrealtype>> right_scale_vec_;
  BlockMatrixType* matrix_;
  // std::shared_ptr<gko::log::BatchConvergence<sunrealtype>> logger_;
  sunrealtype avg_iter_count_;
  sunrealtype sum_of_avg_iters_;
  sunrealtype stddev_iter_count_;
  bool log_res_norm_;
  // std::shared_ptr<const GkoBatchVecType> res_norm_;
  SUNLogger sunlogger_;

  void initSUNLinSol(SUNContext sunctx)
  {
    // Attach function pointers for SUNLinearSolver
    using this_type = BlockLinearSolver<GkoBatchSolverType, BlockMatrixType>;

    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    sunlinsol_->sunctx  = sunctx;

    std::memset(sunlinsol_->ops, 0, sizeof(_generic_SUNLinearSolver_Ops));
    sunlinsol_->ops->gettype           = SUNLinSolGetType_GinkgoBlock;
    sunlinsol_->ops->getid             = SUNLinSolGetID_GinkgoBlock;
    sunlinsol_->ops->setscalingvectors = SUNLinSolSetScalingVectors_GinkgoBlock<this_type>;
    sunlinsol_->ops->initialize        = SUNLinSolInitialize_GinkgoBlock<this_type>;
    sunlinsol_->ops->setup             = SUNLinSolSetup_GinkgoBlock<this_type, BlockMatrixType>;
    sunlinsol_->ops->solve             = SUNLinSolSolve_GinkgoBlock<this_type>;
    sunlinsol_->ops->free              = SUNLinSolFree_GinkgoBlock;
    sunlinsol_->ops->numiters          = SUNLinSolNumIters_GinkgoBlock<this_type>;
  }

  static constexpr int default_max_iters_ = 500;
  static constexpr int default_restart_   = 0;
};

}
} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGOBLOCK_HPP
