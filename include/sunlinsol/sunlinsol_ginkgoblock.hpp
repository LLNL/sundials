/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------*/

#include <cmath>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/priv/sundials_logger_impl.h>
#include <sundials/sundials_core.hpp>

#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

#if (GKO_VERSION_MAJOR < 1) || (GKO_VERSION_MAJOR == 1 && GKO_VERSION_MINOR < 9)
#error "Ginkgo 1.9.0 or later is required."
#endif

#ifndef _SUNLINSOL_GINKGOBLOCK_HPP
#define _SUNLINSOL_GINKGOBLOCK_HPP

namespace sundials {
namespace ginkgo {

template<class GkoBatchSolverType, class GkoBatchMatType>
class BlockLinearSolver;

namespace impl {

inline SUNLinearSolver_Type SUNLinSolGetType_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

inline SUNLinearSolver_ID SUNLinSolGetID_GinkgoBlock(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGOBLOCK;
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolInitialize_GinkgoBlock(SUNLinearSolver S)
{
  return SUN_SUCCESS;
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S,
                                                  N_Vector s1, N_Vector s2)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  SUNFunctionBegin(S->sunctx);
  SUNAssert(s1 || s2, SUN_ERR_ARG_INCOMPATIBLE);
  solver->SetScalingVectors(s1, s2);
  return SUN_SUCCESS;
}

template<class GkoBatchLinearSolverType, class GkoBatchMatType>
int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->Setup(static_cast<BlockMatrix<GkoBatchMatType>*>(A->content));
}

template<class GkoBatchLinearSolverType>
int SUNLinSolSolve_GinkgoBlock(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->Solve(b, x, tol);
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolFree_GinkgoBlock(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  delete solver; // NOLINT(cppcoreguidelines-owning-memory)
  return SUN_SUCCESS;
}

template<class GkoBatchLinearSolverType>
int SUNLinSolNumIters_GinkgoBlock(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return static_cast<int>(std::round(solver->SumAvgNumIters()));
}

inline gko::array<sunrealtype> WrapBatchScalingArray(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_blocks,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  auto xvec_len{N_VGetLength(x)};
  return gko::array<sunrealtype>(gko_exec,
                                 gko::array<sunrealtype>::view(gko_exec,
                                                               xvec_len, x_arr));
}

inline std::unique_ptr<GkoBatchVecType> WrapAsMultiVector(
  gko::array<sunrealtype>& arr, gko::size_type num_blocks)
{
  return GkoBatchVecType::create(arr.get_executor(),
                                 gko::batch_dim<2>(num_blocks,
                                                   gko::dim<2>(arr.get_size() /
                                                                 num_blocks,
                                                               1)),
                                 arr);
}

} // namespace impl

template<class GkoBatchSolverType, class GkoBatchMatType>
class BlockLinearSolver : public sundials::impl::BaseLinearSolver,
                          public ConvertibleTo<SUNLinearSolver>
{
public:
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    int max_iters, gko::size_type num_blocks, SUNContext sunctx)
    : sundials::impl::BaseLinearSolver(sunctx),
      gko_exec_(std::move(gko_exec)),
      tolerance_type_(tolerance_type),
      precon_factory_(std::move(precon_factory)),
      solver_factory_(nullptr),
      solver_(nullptr),
      col_scale_vec_(gko_exec_),
      row_scale_vec_(gko_exec_),
      matrix_(nullptr),
      logger_(nullptr),
      num_blocks_(num_blocks),
      max_iters_(max_iters),
      res_norm_array_(gko_exec_->get_master(), num_blocks_),
      iter_count_array_(gko_exec_->get_master(), num_blocks_),
      avg_iter_count_(sunrealtype{0.0}),
      sum_of_avg_iters_(sunrealtype{0.0}),
      stddev_iter_count_(sunrealtype{0.0}),
      previous_max_res_norm_(sunrealtype{0.0}),
      ones_(nullptr),
      s2inv_(nullptr),
      scaling_(false)
  {
    initSUNLinSol(sunctx);
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, gko::batch::stop::tolerance_type::absolute,
                        nullptr, default_max_iters_, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, tolerance_type, nullptr, default_max_iters_,
                        num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, gko::batch::stop::tolerance_type::absolute,
                        precon_factory, default_max_iters_, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    int max_iters, gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, gko::batch::stop::tolerance_type::absolute,
                        nullptr, max_iters, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    int max_iters, gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, tolerance_type, nullptr, max_iters,
                        num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    int max_iters, gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, gko::batch::stop::tolerance_type::absolute,
                        precon_factory, max_iters, num_blocks, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BlockLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    gko::size_type num_blocks, SUNContext sunctx)
    : BlockLinearSolver(gko_exec, tolerance_type, precon_factory,
                        default_max_iters_, num_blocks, sunctx)
  {}

  /// Get the ``gko::Executor`` associated with the Ginkgo matrix
  std::shared_ptr<const gko::Executor> GkoExec() const { return gko_exec_; }

  // Override the ConvertibleTo methods

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() override { return object_.get(); }

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() const override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver Convert() override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver Convert() const override { return object_.get(); }

  /// Get the underlying Ginkgo solver
  /// \note This will be `nullptr` until the linear solver setup phase.
  GkoBatchSolverType* GkoSolver() { return solver_.get(); }

  /// Average number of iterations across the batches during the last solve.
  sunrealtype AvgNumIters() const { return avg_iter_count_; }

  /// Standard deviation of the number of iterations across the batches during the last solve.
  sunrealtype StddevNumIters() const { return stddev_iter_count_; }

  /// Running sum of the average number of iterations in this solvers lifetime.
  sunrealtype SumAvgNumIters() const { return sum_of_avg_iters_; }

  /// Sets the left and right scaling vectors to be used.
  void SetScalingVectors(N_Vector s1, N_Vector s2)
  {
    if (s1 && s2)
    {
      col_scale_vec_ =
        std::move(impl::WrapBatchScalingArray(GkoExec(), num_blocks_, s1));

      if (!ones_.Convert())
      {
        ones_ = sundials::experimental::NVectorView(N_VClone(s2));
        N_VConst(sunrealtype{1.0}, ones_);
      }

      if (!s2inv_.Convert())
      {
        s2inv_ = sundials::experimental::NVectorView(N_VClone(s2));
      }

      N_VDiv(ones_, s2,
             s2inv_); // Need to provide s2inv to match the SUNLinearSolver API

      row_scale_vec_ =
        std::move(impl::WrapBatchScalingArray(GkoExec(), num_blocks_, s2inv_));

      scaling_ = true;
    }
    else { scaling_ = false; }
  }

  int Setup(BlockMatrix<GkoBatchMatType>* A)
  {
    if (num_blocks_ != A->NumBlocks()) { return SUN_ERR_ARG_OUTOFRANGE; }

    matrix_ = A;

    if (!scaling_)
    {
      solver_factory_ = GkoBatchSolverType::build()             //
                          .with_max_iterations(max_iters_)      //
                          .with_tolerance(SUN_UNIT_ROUNDOFF)    //
                          .with_tolerance_type(tolerance_type_) //
                          .with_preconditioner(precon_factory_) //
                          .on(GkoExec());

      solver_ = solver_factory_->generate(matrix_->GkoMtx());
    }

    return SUN_SUCCESS;
  }

  int Solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    SUNLogInfo(sunLogger(), "linear-solver", "solver = ginkgo (batched)");

    if (scaling_)
    {
      solver_factory_ = GkoBatchSolverType::build()             //
                          .with_max_iterations(max_iters_)      //
                          .with_tolerance(tol)                  //
                          .with_tolerance_type(tolerance_type_) //
                          .with_preconditioner(precon_factory_) //
                          .on(GkoExec());

      // \tilde{A} = S_1 A S_2^{-1}
      matrix_->GkoMtx()->scale(row_scale_vec_, col_scale_vec_);

      solver_ = solver_factory_->generate(matrix_->GkoMtx());
    }
    else { solver_->reset_tolerance(tol); }

    if (!logger_)
    {
      logger_ =
        gko::share(gko::batch::log::BatchConvergence<sunrealtype>::create());
    }
    solver_->add_logger(logger_);

    std::unique_ptr<GkoBatchVecType> x_vec{
      impl::WrapBatchVector(GkoExec(), num_blocks_, x)};
    std::unique_ptr<GkoBatchVecType> b_vec{
      impl::WrapBatchVector(GkoExec(), num_blocks_, b)};

    if (scaling_)
    {
      // \tilde{b} = S_1 b
      b_vec->scale(impl::WrapAsMultiVector(col_scale_vec_, num_blocks_));
    }

    // \tilde{x} = \tilde{A}^{-1} \tilde{b}
    [[maybe_unused]] gko::batch::BatchLinOp* result =
      solver_->apply(b_vec.get(), x_vec.get());

    if (scaling_)
    {
      // x = S_2^{-1} \tilde{x}
      x_vec->scale(impl::WrapAsMultiVector(row_scale_vec_, num_blocks_));
    }

    // Check if any batch entry did not reach the tolerance.
    bool at_least_one_did_not_converge{false};
    bool max_res_norm_did_not_reduce{false};
    sunrealtype max_res_norm{0.0};
    sunrealtype min_res_norm{SUN_BIG_REAL};
    gko_exec_->get_master()->copy_from(gko_exec_, num_blocks_,
                                       logger_->get_residual_norm().get_const_data(),
                                       res_norm_array_.get_data());
    auto res_norm = res_norm_array_.get_data();

    for (gko::size_type i = 0; i < num_blocks_; i++)
    {
      max_res_norm =
        std::max(max_res_norm,
                 res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_res_norm =
        std::min(min_res_norm,
                 res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    if (max_res_norm > tol) { at_least_one_did_not_converge = true; }
    if (max_res_norm >= previous_max_res_norm_)
    {
      max_res_norm_did_not_reduce = true;
    }

    // Compute the average number of iterations across all batch entries
    // as well as the maximum and minimum iteration counts.
    gko_exec_->get_master()->copy_from(gko_exec_, num_blocks_,
                                       logger_->get_num_iterations().get_const_data(),
                                       iter_count_array_.get_data());
    auto iter_count = iter_count_array_.get_data();

    sunindextype max_iter_count{0};
    sunindextype min_iter_count{max_iters_};
    avg_iter_count_ = 0.0;
    for (gko::size_type i = 0; i < num_blocks_; i++)
    {
      avg_iter_count_ +=
        iter_count[i]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      max_iter_count =
        std::max(max_iter_count,
                 iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_iter_count =
        std::min(min_iter_count,
                 iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    avg_iter_count_ /= static_cast<sunrealtype>(num_blocks_);
    sum_of_avg_iters_ += avg_iter_count_;

    // Compute the std. dev. in iteration count across all batch entries.
    // This helps us understand how varied (in difficulty to solve) the entries are.
    stddev_iter_count_ = 0.0;
    for (gko::size_type i = 0; i < num_blocks_; i++)
    {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      stddev_iter_count_ += std::pow(std::abs(iter_count[i] - avg_iter_count_),
                                     2);
    }
    stddev_iter_count_ =
      std::sqrt(stddev_iter_count_ / static_cast<sunrealtype>(num_blocks_));

    int retval{0};
    if (at_least_one_did_not_converge && max_res_norm_did_not_reduce)
    {
      retval = SUNLS_CONV_FAIL;
    }
    else if (at_least_one_did_not_converge) { retval = SUNLS_RES_REDUCED; }
    else { retval = SUN_SUCCESS; }

    SUNLogInfo(sunLogger(), "linear-solver",
               "avg. iter count = %.16g, stddev. iter count = %.16g, "
               "max iter count = %d, "
               "min iter count = %d, max res. norm = %.16g, min res. "
               "norm = %.16g, tol = %.16g,"
               " return code = %d",
               avg_iter_count_, stddev_iter_count_, max_iter_count,
               min_iter_count, max_res_norm, min_res_norm, tol, retval);

    return retval;
  }

private:
  std::shared_ptr<const gko::Executor> gko_exec_;
  gko::batch::stop::tolerance_type tolerance_type_;
  std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory_;
  std::unique_ptr<typename GkoBatchSolverType::Factory> solver_factory_;
  std::unique_ptr<GkoBatchSolverType> solver_;
  gko::array<sunrealtype> col_scale_vec_;
  gko::array<sunrealtype> row_scale_vec_;
  BlockMatrix<GkoBatchMatType>* matrix_;
  std::shared_ptr<gko::batch::log::BatchConvergence<sunrealtype>> logger_;
  gko::size_type num_blocks_;
  int max_iters_;
  gko::array<sunrealtype> res_norm_array_;
  gko::array<sunindextype> iter_count_array_;
  sunrealtype avg_iter_count_;
  sunrealtype sum_of_avg_iters_;
  sunrealtype stddev_iter_count_;
  sunrealtype previous_max_res_norm_;
  sundials::experimental::NVectorView ones_, s2inv_;
  bool scaling_;

  void initSUNLinSol(SUNContext sunctx)
  {
    // Attach function pointers for SUNLinearSolver
    using this_type = BlockLinearSolver<GkoBatchSolverType, GkoBatchMatType>;

    object_->content = this;
    object_->ops     = object_ops_.get();
    object_->sunctx  = sunctx;

    object_->ops->gettype = impl::SUNLinSolGetType_GinkgoBlock;
    object_->ops->getid   = impl::SUNLinSolGetID_GinkgoBlock;
    object_->ops->setscalingvectors =
      impl::SUNLinSolSetScalingVectors_GinkgoBlock<this_type>;
    object_->ops->initialize = impl::SUNLinSolInitialize_GinkgoBlock<this_type>;
    object_->ops->setup =
      impl::SUNLinSolSetup_GinkgoBlock<this_type, GkoBatchMatType>;
    object_->ops->solve    = impl::SUNLinSolSolve_GinkgoBlock<this_type>;
    object_->ops->free     = impl::SUNLinSolFree_GinkgoBlock<this_type>;
    object_->ops->numiters = impl::SUNLinSolNumIters_GinkgoBlock<this_type>;
  }

  SUNLogger sunLogger()
  {
    SUNLogger log{nullptr};
    SUNContext_GetLogger(object_->sunctx, &log);
    return log;
  }

  static constexpr int default_max_iters_ = 500;
  static constexpr int default_restart_   = 0;
};

} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGOBLOCK_HPP
