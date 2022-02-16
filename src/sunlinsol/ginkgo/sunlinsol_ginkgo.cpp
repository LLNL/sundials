#include <sunlinsol/sunlinsol_ginkgo.hpp>

using namespace sundials;

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;

#define GET_CONTENT(S) \
  ((ginkgo::LinearSolver<gko::solver::Cg<sunrealtype>, GkoDenseMat>*)S->content)
