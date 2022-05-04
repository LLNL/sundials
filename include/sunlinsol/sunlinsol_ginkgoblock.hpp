
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
  return SUNMAT_SUCCESS;
}

template<typename GkoBatchLinearSolverType>
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_GinkgoBlock(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  // auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return -1;
}

template<typename GkoBatchLinearSolverType, typename BlockMatrixType>
SUNDIALS_EXPORT int SUNLinSolSetup_GinkgoBlock(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  solver->setup(static_cast<BlockMatrixType*>(A->content));
  return SUNLS_SUCCESS;
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

class BlockLinearSolverInterface {
public:
  virtual gko::BatchLinOp* solve(N_Vector b, N_Vector x, sunrealtype tol) = 0;
  virtual std::shared_ptr<const gko::Executor> gkoexec() const       = 0;
  virtual SUNLinearSolver get()                                      = 0;
  virtual SUNLinearSolver get() const                                = 0;
  virtual operator SUNLinearSolver()                                 = 0;
  virtual operator SUNLinearSolver() const                           = 0;
};

template<class GkoBatchSolverType, class BlockMatrixType>
class BlockLinearSolver : public BlockLinearSolverInterface {
public:
  BlockLinearSolver(std::shared_ptr<typename GkoBatchSolverType::Factory> gko_solver_factory, SUNContext sunctx)
      : gko_solver_factory_(gko_solver_factory), sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
        sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>())
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

  std::shared_ptr<const gko::Executor> gkoexec() const { return gko_solver_->get_executor(); }

  std::shared_ptr<typename GkoBatchSolverType::Factory> gkofactory() { return gko_solver_factory_; }

  GkoBatchSolverType* gkosolver() { return gko_solver_.get(); }

  operator SUNLinearSolver() { return sunlinsol_.get(); }

  operator SUNLinearSolver() const { return sunlinsol_.get(); }

  SUNLinearSolver get() { return sunlinsol_.get(); }

  SUNLinearSolver get() const { return sunlinsol_.get(); }

  GkoBatchSolverType* setup(BlockMatrixType* A)
  {
    gko_solver_ = gko_solver_factory_->generate(A->gkomtx());
    return gko_solver_.get();
  }

  gko::BatchLinOp* solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    gko::BatchLinOp* result = nullptr;
    sunindextype num_blocks = gkosolver()->get_num_batch_entries();

    if (x != b) {
      auto x_vec = WrapBatchVector(gkoexec(), num_blocks, x);
      auto b_vec = WrapBatchVector(gkoexec(), num_blocks, b);

      // x = A^{-1} b
      result = gkosolver()->apply(b_vec.get(), x_vec.get());
    } else {
      auto x_vec = WrapBatchVector(gkoexec(), num_blocks, x);

      // x = A^{-1} x
      result = gkosolver()->apply(x_vec.get(), x_vec.get());
    }

    return result;
  }

private:
  std::shared_ptr<typename GkoBatchSolverType::Factory> gko_solver_factory_;
  std::unique_ptr<GkoBatchSolverType> gko_solver_;
  std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<struct _generic_SUNLinearSolver_Ops> sunlinsol_ops_;
};

} // namespace ginkgo
} // namespace sundials

#endif
