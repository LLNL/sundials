#include <sunlinsol/sunlinsol_ginkgo.hpp>

using namespace sundials;

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;

#define GET_CONTENT(S) ((ginkgo::LinearSolver<gko::solver::Cg<>, GkoDenseMat> *) S->content)

SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S)
{
   return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S)
{
   return SUNLINEARSOLVER_GINKGO;
}

int SUNLinSolInitialize_Ginkgo(SUNLinearSolver S)
{
   return SUNMAT_SUCCESS;
}

int SUNLinSolSetATimes_Ginkgo(SUNLinearSolver S, void* A_data,
                              SUNATimesFn ATimes);
int SUNLinSolSetPreconditioner_Ginkgo(SUNLinearSolver S,
                                      void* P_data,
                                      SUNPSetupFn Pset,
                                      SUNPSolveFn Psol);
int SUNLinSolSetScalingVectors_Ginkgo(SUNLinearSolver S,
                                      N_Vector s1,
                                      N_Vector s2);
int SUNLinSolSetZeroGuess_Ginkgo(SUNLinearSolver S,
                                 booleantype onff);


int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
   auto solver = GET_CONTENT(S);
   solver->generate(static_cast<ginkgo::Matrix<GkoDenseMat>*>(A->content)->gkomtx());
   return SUNMAT_SUCCESS;
}

int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A,
                          N_Vector x, N_Vector b, realtype tol)
{

}

int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S);
realtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S);
N_Vector SUNLinSolResid_Ginkgo(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_Ginkgo(SUNLinearSolver S);
int SUNLinSolSpace_Ginkgo(SUNLinearSolver S,
                          long int *lenrwLS,
                          long int *leniwLS);
int SUNLinSolFree_Ginkgo(SUNLinearSolver S);