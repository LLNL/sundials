#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <memory>

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_linearsolver.hpp>
#include <sundials/sundials_matrix.hpp>
#include <sundials/sundials_nvector.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "ida/ida_impl.h"

int res(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector residual,
        void* user_data)
{
  NV_Ith_S(residual, 0) = NV_Ith_S(ydot, 0) - std::cos(t);
  return 0;
}

TEST(IDAGetInterpData, ReturnsInterpolationDataForLastStep)
{
  using namespace sundials::experimental;

  SUNErrCode sunerr;

  const sunrealtype t0    = 0.0;
  const sunrealtype y0    = 0.0;
  const sunrealtype ydot0 = 1.0;

  sundials::Context sunctx(SUN_COMM_NULL);

  NVectorView y            = N_VNew_Serial(1, sunctx);
  NV_Ith_S(y.Convert(), 0) = y0;

  NVectorView ydot            = N_VNew_Serial(1, sunctx);
  NV_Ith_S(ydot.Convert(), 0) = ydot0;

  std::unique_ptr<void, std::function<void(void*)>> ida_mem(IDACreate(sunctx),
                                                            [](void* p)
                                                            { IDAFree(&p); });

  ASSERT_EQ(IDAInit(ida_mem.get(), res, t0, y, ydot), 0);
  ASSERT_EQ(IDASStolerances(ida_mem.get(), 1e-6, 1e-10), 0);

  SUNMatrixView A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolverView LS = SUNLinSol_Dense(y, A, sunctx);
  ASSERT_EQ(IDASetLinearSolver(ida_mem.get(), LS, A), 0);

  const int N_STEPS      = 10;
  const sunrealtype tout = 1.0;
  sunrealtype t;
  for (int n_step = 0; n_step < N_STEPS; ++n_step)
  {
    ASSERT_EQ(IDASolve(ida_mem.get(), tout, &t, y, ydot, IDA_ONE_STEP), 0);

    N_Vector* phi;
    sunrealtype* psi;
    int kused;
    sunrealtype hused;
    sunrealtype tn;
    ASSERT_EQ(IDAGetInterpData(ida_mem.get(), &phi, &psi, &kused, &hused, &tn),
              0);

    int ida_kused         = -1;
    sunrealtype ida_hused = -1.0;
    sunrealtype ida_tn    = -1.0;
    IDAGetLastOrder(ida_mem.get(), &ida_kused);
    IDAGetLastStep(ida_mem.get(), &ida_hused);
    IDAGetCurrentTime(ida_mem.get(), &ida_tn);

    auto ida_phi = static_cast<IDAMemRec*>(ida_mem.get())->ida_phi;
    auto ida_psi = static_cast<IDAMemRec*>(ida_mem.get())->ida_psi;

    ASSERT_EQ(phi, ida_phi);
    ASSERT_EQ(psi, ida_psi);
    ASSERT_EQ(kused, ida_kused);
    ASSERT_EQ(hused, ida_hused);
    ASSERT_EQ(tn, ida_tn);
  }
}
