#include <cmath>
#include <memory>
#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_logger.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>
#include <utility>

#ifndef _SUNLINSOL_GINKGOBLOCK_HPP
#define _SUNLINSOL_GINKGOBLOCK_HPP

namespace sundials {
namespace ginkgo {

template<class GkoBatchSolverType, class BlockMatrixType>
class BlockLinearSolver;

SUNLinearSolver_Type SUNLinSolGetType_GinkgoBlock(SUNLinearSolver S) { return SUNLINEARSOLVER_MATRIX_ITERATIVE; }

SUNLinearSolver_ID SUNLinSolGetID_GinkgoBlock(SUNLinearSolver S) { return SUNLINEARSOLVER_GINKGOBLOCK; }

template<typename GkoBatchLinearSolverType>
int SUNLinSolInitialize_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLS_SUCCESS;
}

template<typename GkoBatchLinearSolverType>
int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  if (!s1 || !s2) {
    return SUNLS_ILL_INPUT;
  }
  solver->setScalingVectors(s1, s2);
  return SUNLS_SUCCESS;
}

template<typename GkoBatchLinearSolverType, typename BlockMatrixType>
int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->setup(static_cast<BlockMatrixType*>(A->content));
}

template<typename GkoBatchLinearSolverType>
int SUNLinSolSolve_GinkgoBlock(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->solve(b, x, tol);
}

template<typename GkoBatchLinearSolverType>
int SUNLinSolFree_GinkgoBlock(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  delete solver; // NOLINT(cppcoreguidelines-owning-memory)
  return SUNLS_SUCCESS;
}

template<typename GkoBatchLinearSolverType>
int SUNLinSolNumIters_GinkgoBlock(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return std::round(solver->sumAvgNumIters());
}

inline std::unique_ptr<gko::matrix::BatchDiagonal<sunrealtype>> WrapBatchDiagMatrix(
    std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks, N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x) : N_VGetArrayPointer(x)};
  auto xvec_len{N_VGetLength(x)};
  auto batch_xvec_size{gko::batch_dim<2>(num_blocks, gko::dim<2>(xvec_len / num_blocks, xvec_len / num_blocks))};
  auto xvec_view{gko::Array<sunrealtype>::view(gko_exec, xvec_len, x_arr)};
  return gko::matrix::BatchDiagonal<sunrealtype>::create(gko_exec, batch_xvec_size, std::move(xvec_view));
}

template<class GkoBatchSolverType, class BlockMatrixType>
class BlockLinearSolver : public ConvertibleTo<SUNLinearSolver> {
public:
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    std::shared_ptr<gko::BatchLinOpFactory> precon_factory, int max_iters, sunindextype num_blocks,
                    SUNContext sunctx)
      : gko_exec_(std::move(gko_exec)), tolerance_type_(tolerance_type), precon_factory_(std::move(precon_factory)),
        solver_factory_(nullptr), solver_(nullptr), sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
        sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>()), left_scale_vec_(nullptr),
        right_scale_vec_(nullptr), matrix_(nullptr), logger_(nullptr), num_blocks_(num_blocks), max_iters_(max_iters)
  {
    initSUNLinSol(sunctx);
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, nullptr, default_max_iters_, num_blocks,
                          sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, nullptr, default_max_iters_, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    std::shared_ptr<gko::BatchLinOpFactory> precon_factory, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, precon_factory, default_max_iters_,
                          num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, int max_iters, sunindextype num_blocks,
                    SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, nullptr, max_iters, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    int max_iters, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, nullptr, max_iters, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, std::shared_ptr<gko::BatchLinOpFactory> precon_factory,
                    int max_iters, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, gko::stop::batch::ToleranceType::absolute, precon_factory, max_iters, num_blocks,
                          sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec, gko::stop::batch::ToleranceType tolerance_type,
                    std::shared_ptr<gko::BatchLinOpFactory> precon_factory, sunindextype num_blocks, SUNContext sunctx)
      : BlockLinearSolver(gko_exec, tolerance_type, precon_factory, default_max_iters_, num_blocks, sunctx)
  {}

  std::shared_ptr<const gko::Executor> gkoExec() const { return gko_exec_; }

  void setEnableScaling(bool onoff)
  {
    if (!onoff) {
      sunlinsol_->ops->setscalingvectors = nullptr;
      left_scale_vec_                    = nullptr;
      right_scale_vec_                   = nullptr;
    } else {
      sunlinsol_->ops->setscalingvectors =
          SUNLinSolSetScalingVectors_GinkgoBlock<BlockLinearSolver<GkoBatchSolverType, BlockMatrixType>>;
    }
  }

  SUNLinearSolver get() override { return sunlinsol_.get(); }
  SUNLinearSolver get() const override { return sunlinsol_.get(); }
  operator SUNLinearSolver() override { return sunlinsol_.get(); }
  operator SUNLinearSolver() const override { return sunlinsol_.get(); }

  sunrealtype avgNumIters() const { return avg_iter_count_; }

  sunrealtype stddevNumIters() const { return stddev_iter_count_; }

  sunrealtype sumAvgNumIters() const { return sum_of_avg_iters_; }

  void setScalingVectors(N_Vector s1, N_Vector s2)
  {
    if (!ones_) {
      ones_ = N_VClone(s2);
      N_VConst(sunrealtype{1.0}, ones_);
    }
    if (!s2inv_) {
      s2inv_ = N_VClone(s2);
    }
    N_VDiv(ones_, s2, s2inv_); // Need to provide s2inv to match the SUNLinearSolver API
    left_scale_vec_  = gko::share(WrapBatchDiagMatrix(gkoExec(), num_blocks_, s1));
    right_scale_vec_ = gko::share(WrapBatchDiagMatrix(gkoExec(), num_blocks_, s2inv_));
  }

  int setup(BlockMatrixType* A)
  {
    if (num_blocks_ != A->numBlocks()) {
      return SUNLS_ILL_INPUT;
    }
    setup_called_ = true;
    matrix_       = A;
    return SUNLS_SUCCESS;
  }

  int solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(sunLogger(), SUN_LOGLEVEL_INFO, "sundials::ginkgo::BlockLinearSolver::solve", "start",
                       "num_blocks = %d, tol = %.16g", num_blocks_, tol);
#endif

    // Ideally we would only build and generate the solver when the matrix
    // changes, but we cannot since the tolerance value is not known until now.
    // TOOD(CJB): determine how expensive this build+generate is and ask Ginkgo for a way to avoid it if it is expensive
    SUNDIALS_MARK_BEGIN(sunProfiler(), "build solver factory");
    solver_factory_ = GkoBatchSolverType::build()                  //
                          .with_max_iterations(max_iters_)         //
                          .with_residual_tol(tol)                  //
                          .with_tolerance_type(tolerance_type_)    //
                          .with_preconditioner(precon_factory_)    //
                          .with_left_scaling_op(left_scale_vec_)   //
                          .with_right_scaling_op(right_scale_vec_) //
                          .on(gkoExec());
    SUNDIALS_MARK_END(sunProfiler(), "build solver factory");

    SUNDIALS_MARK_BEGIN(sunProfiler(), "generate solver");
    solver_ = solver_factory_->generate(matrix_->gkomtx());
    SUNDIALS_MARK_END(sunProfiler(), "generate solver");

    SUNDIALS_MARK_BEGIN(sunProfiler(), "add logger");
    if (!logger_) {
      logger_ = gko::share(gko::log::BatchConvergence<sunrealtype>::create(gkoExec()));
    }
    solver_->add_logger(logger_);
    SUNDIALS_MARK_END(sunProfiler(), "add logger");

    gko::BatchLinOp* result{nullptr};
    if (x != b) {
      SUNDIALS_MARK_BEGIN(sunProfiler(), "Wrap vector(s) for solve");
      std::unique_ptr<GkoBatchVecType> x_vec = WrapBatchVector(gkoExec(), num_blocks_, x);
      std::unique_ptr<GkoBatchVecType> b_vec = WrapBatchVector(gkoExec(), num_blocks_, b);
      SUNDIALS_MARK_END(sunProfiler(), "Wrap vector(s) for solve");

      // x = A'^{-1} diag(left) b
      SUNDIALS_MARK_BEGIN(sunProfiler(), "solver apply");
      // printf("b:\n"); N_VPrint(b);
      result = solver_->apply(b_vec.get(), x_vec.get());
      // printf("x:\n"); N_VPrint(x);
      SUNDIALS_MARK_END(sunProfiler(), "solver apply");
    } else {
      SUNDIALS_MARK_BEGIN(sunProfiler(), "Wrap vector(s) for solve");
      std::unique_ptr<GkoBatchVecType> x_vec = WrapBatchVector(gkoExec(), num_blocks_, x);
      SUNDIALS_MARK_END(sunProfiler(), "Wrap vector(s) for solve");

      // x = A^'{-1} diag(right) x
      SUNDIALS_MARK_BEGIN(sunProfiler(), "solver apply");
      // printf("b:\n"); N_VPrint(b);
      result = solver_->apply(x_vec.get(), x_vec.get());
      // printf("x:\n"); N_VPrint(x);
      SUNDIALS_MARK_END(sunProfiler(), "solver apply");
    }

    SUNDIALS_MARK_BEGIN(sunProfiler(), "check residual norm");
    // Check if any batch entry did not reach the tolerance.
    bool at_least_one_did_not_converge{false};
    bool max_res_norm_did_not_reduce{false};
    sunrealtype max_res_norm{0.0};
    sunrealtype min_res_norm{SUN_BIG_REAL};
    auto res_norm{logger_->get_residual_norm()->get_const_values()};
    for (int i = 0; i < num_blocks_; i++) {
      max_res_norm = std::max(max_res_norm, res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_res_norm = std::min(min_res_norm, res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    if (max_res_norm > tol) {
      at_least_one_did_not_converge = true;
    }
    if (max_res_norm >= previous_max_res_norm_) {
      max_res_norm_did_not_reduce = true;
    }
    SUNDIALS_MARK_END(sunProfiler(), "check residual norm");

    SUNDIALS_MARK_BEGIN(sunProfiler(), "check num iters");
    // Compute the average number of iterations across all batch entries
    // as well as the maximum and minimum iteration counts.
    auto iter_count{logger_->get_num_iterations().get_const_data()};
    int max_iter_count{0};
    int min_iter_count{max_iters_};
    avg_iter_count_ = 0.0;
    for (int i = 0; i < num_blocks_; i++) {
      avg_iter_count_ += iter_count[i]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      max_iter_count = std::max(max_iter_count, iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_iter_count = std::min(min_iter_count, iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    avg_iter_count_ /= num_blocks_;
    sum_of_avg_iters_ += avg_iter_count_;

    // Compute the std. dev. in iteration count across all batch entries.
    // This helps us understand how varied (in difficulty to solve) the entries are.
    stddev_iter_count_ = 0.0;
    for (int i = 0; i < num_blocks_; i++) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      stddev_iter_count_ += std::pow(std::abs(iter_count[i] - avg_iter_count_), 2);
    }
    stddev_iter_count_ = std::sqrt(stddev_iter_count_ / num_blocks_);
    SUNDIALS_MARK_END(sunProfiler(), "check num iters");

    int retval{0};
    if (at_least_one_did_not_converge && max_res_norm_did_not_reduce) {
      // printf("return SUNLS_CONV_FAIL\n");
      retval = SUNLS_CONV_FAIL;
    } else if (at_least_one_did_not_converge) {
      // printf("retval = SUNLS_RES_REDUCED\n");
      retval = SUNLS_RES_REDUCED;
    } else {
      // printf("retval = SUNLS_SUCCESS\n");
      retval = SUNLS_SUCCESS;
    }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    SUNLogger_QueueMsg(sunLogger(), SUN_LOGLEVEL_INFO, "sundials::ginkgo::BlockLinearSolver::solve", "end-solve",
                       "avg. iter count = %.16g, stddev. iter count = %.16g, max iter count = %d, "
                       "min iter count = %d, max res. norm = %.16g, min res. norm = %.16g, tol = %.16g,"
                       " return code = %d",
                       avg_iter_count_, stddev_iter_count_, max_iter_count, min_iter_count, max_res_norm, min_res_norm,
                       tol, retval);
#endif

    return retval;
  }

private:
  std::shared_ptr<const gko::Executor> gko_exec_;
  gko::stop::batch::ToleranceType tolerance_type_;
  std::shared_ptr<gko::BatchLinOpFactory> precon_factory_;
  std::unique_ptr<typename GkoBatchSolverType::Factory> solver_factory_;
  std::unique_ptr<GkoBatchSolverType> solver_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<struct _generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  std::shared_ptr<gko::matrix::BatchDiagonal<sunrealtype>> left_scale_vec_;
  std::shared_ptr<gko::matrix::BatchDiagonal<sunrealtype>> right_scale_vec_;
  BlockMatrixType* matrix_;
  std::shared_ptr<gko::log::BatchConvergence<sunrealtype>> logger_;
  sunindextype num_blocks_{};
  int max_iters_{};
  bool setup_called_{};
  sunrealtype avg_iter_count_{};
  sunrealtype sum_of_avg_iters_{};
  sunrealtype stddev_iter_count_{};
  sunrealtype previous_max_res_norm_{};
  N_Vector ones_{}, s2inv_{};

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
    sunlinsol_->ops->free              = SUNLinSolFree_GinkgoBlock<this_type>;
    sunlinsol_->ops->numiters          = SUNLinSolNumIters_GinkgoBlock<this_type>;
  }

  SUNLogger sunLogger()
  {
    SUNLogger log{nullptr};
    SUNContext_GetLogger(sunlinsol_->sunctx, &log);
    return log;
  }

  SUNProfiler sunProfiler()
  {
    SUNProfiler prof{nullptr};
    SUNContext_GetProfiler(sunlinsol_->sunctx, &prof);
    return prof;
  }

  static constexpr int default_max_iters_ = 500;
  static constexpr int default_restart_   = 0;
};

} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGOBLOCK_HPP
