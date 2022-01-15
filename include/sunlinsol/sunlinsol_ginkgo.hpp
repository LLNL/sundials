
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

template<typename GkoSolverType, typename GkoMatType>
class LinearSolver
{
public:
   LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory, SUNContext sunctx)
      : gkosolver_factory_(gko_solver_factory),
        sunlinsol_(SUNLinSolNewEmpty(sunctx))
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

   std::shared_ptr<GkoSolverType> generate(std::shared_ptr<GkoMatType> A)
   {
      return gko::share(gkosolver_factory_->generate(A));
   }

private:
   std::shared_ptr<typename GkoSolverType::Factory> gkosolver_factory_;
   std::unique_ptr<struct _generic_SUNLinearSolver> sunlinsol_;

};


}// namespace ginkgo
}// namespace sundials

#endif
