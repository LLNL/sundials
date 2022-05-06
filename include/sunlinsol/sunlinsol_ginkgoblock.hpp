
#include <cstring>
#include <memory>
#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgoblock.hpp>

#ifndef _SUNLINSOL_GINKGOBLOCK_HPP
#define _SUNLINSOL_GINKGOBLOCK_HPP

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_GinkgoBlock(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_GinkgoBlock(SUNLinearSolver S, void* P_data, SUNPSetupFn Pset,
                                                           SUNPSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S, N_Vector s1, N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_GinkgoBlock(SUNLinearSolver S, booleantype onff);
SUNDIALS_EXPORT int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_GinkgoBlock(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_GinkgoBlock(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_GinkgoBlock(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_GinkgoBlock(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

namespace sundials {
namespace ginkgo {

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
  solver->solve(b, x, tol);
  return SUNLS_SUCCESS;
}

SUNDIALS_EXPORT
int SUNLinSolFree_GinkgoBlock(SUNLinearSolver S) { return SUNLS_SUCCESS; }

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
        left_scale_vec_(nullptr), right_scale_vec_(nullptr)
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

  void setScalingVectors(N_Vector s1, N_Vector s2)
  {
    left_scale_vec_  = WrapBatchDiagMatrix(gkoExec(), num_blocks_, s1);
    right_scale_vec_ = WrapBatchDiagMatrix(gkoExec(), num_blocks_, s2);
  }

  int setup(BlockMatrixType* A)
  {
    if (num_blocks_ != A->numBlocks()) {
      return SUNLS_ILL_INPUT;
    }
    matrix_ = std::shared_ptr<BlockMatrixType>(A);
    return SUNLS_SUCCESS;
  }

  gko::BatchLinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    gko::BatchLinOp* result = nullptr;

    auto solver_factory = GkoBatchSolverType::build()                                        //
                              .with_max_iterations(max_iters_)                               //
                              .with_residual_tol(tol)                                        //
                              .with_tolerance_type(tolerance_type_)                          //
                              .with_preconditioner(precon_factory_)           //
                              .with_left_scaling_op(left_scale_vec_.get())                   //
                              .with_right_scaling_op(right_scale_vec_.get())                 //
                              .on(gkoExec());

    auto solver = solver_factory->generate(matrix_->gkomtx());

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

    return result;
  }

private:
  std::shared_ptr<const gko::Executor> gko_exec_;
  gko::stop::batch::ToleranceType tolerance_type_;
  std::shared_ptr<gko::BatchLinOpFactory> precon_factory_;
  int max_iters_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<struct _generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  sunindextype num_blocks_;
  std::unique_ptr<gko::BatchLinOp> left_scale_vec_;
  std::unique_ptr<gko::BatchLinOp> right_scale_vec_;
  std::shared_ptr<BlockMatrixType> matrix_;

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
  }

  static constexpr int default_max_iters_ = 50;
  static constexpr int default_restart_   = 0;
};

} // namespace ginkgo
} // namespace sundials

#endif
