/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_nvector_base.hpp>

namespace
{
class TestVector final : public sundials::NVectorBase<TestVector>
{
public:
  static int live;
  static int destroyed;
  static int clones;
  static int fail_after;

  explicit TestVector(SUNContext ctx) : NVectorBase(ctx) { ++live; }
  ~TestVector() override
  {
    --live;
    ++destroyed;
  }

  TestVector* Clone(SUNContext ctx) const override
  {
    if (fail_after >= 0 && clones++ >= fail_after) return nullptr;
    return new TestVector(ctx);
  }

  sunindextype GetLength() const override { return 1; }
  void LinearSum(sunrealtype, const TestVector&, sunrealtype, const TestVector&) override {}
  void Const(sunrealtype) override {}
  void Scale(sunrealtype, const TestVector&) override {}
  sunrealtype DotProd(const TestVector&) const override { return 0; }
  sunrealtype MaxNorm() const override { return 0; }

public:
  SUNErrCode LinearCombination(ConstReals, ConstVectors) override
  {
    throw std::runtime_error("bridge test exception");
  }
};

int TestVector::live       = 0;
int TestVector::destroyed  = 0;
int TestVector::clones     = 0;
int TestVector::fail_after = -1;

struct ErrorCapture
{
  int line = 0;
  std::string func;
  std::string message;
  SUNErrCode code = SUN_SUCCESS;
};

void CaptureError(int line, const char* func, const char*, const char* msg,
                  SUNErrCode code, void* data, SUNContext)
{
  auto* capture   = static_cast<ErrorCapture*>(data);
  capture->line   = line;
  capture->func   = func;
  capture->message = msg;
  capture->code   = code;
}

class NVectorBaseTest : public testing::Test
{
protected:
  void SetUp() override
  {
    ASSERT_EQ(SUNContext_Create(SUN_COMM_NULL, &ctx), SUN_SUCCESS);
    SUNContext_ClearErrHandlers(ctx);
    ASSERT_EQ(SUNContext_PushErrHandler(ctx, CaptureError, &error), SUN_SUCCESS);
    TestVector::live = TestVector::destroyed = TestVector::clones = 0;
    TestVector::fail_after = -1;
  }

  void TearDown() override
  {
    EXPECT_EQ(TestVector::live, 0);
    SUNContext_Free(&ctx);
  }

  SUNContext ctx = nullptr;
  ErrorCapture error;
};

TEST_F(NVectorBaseTest, UserHandleAndCloneOwnership)
{
  {
    TestVector user(ctx);
    N_Vector handle = user.GetNVector();
    error.line    = 123;
    error.func    = "before destroy";
    error.message = "before destroy";
    error.code    = SUN_SUCCESS;
    N_VDestroy(handle);
    EXPECT_EQ(SUNContext_GetLastError(ctx), SUN_SUCCESS);
    EXPECT_EQ(error.line, 123);
    EXPECT_EQ(error.func, "before destroy");
    EXPECT_EQ(error.message, "before destroy");
    EXPECT_EQ(error.code, SUN_SUCCESS);
    EXPECT_EQ(TestVector::live, 1);

    N_Vector clone = N_VClone(handle);
    ASSERT_NE(clone, nullptr);
    EXPECT_EQ(TestVector::live, 2);
    N_VDestroy(clone);
    EXPECT_EQ(TestVector::live, 1);
    EXPECT_EQ(TestVector::destroyed, 1);
  }
  EXPECT_EQ(TestVector::destroyed, 2);
}

TEST_F(NVectorBaseTest, CloneArrayFailureCleansEarlierClones)
{
  TestVector user(ctx);
  // The clone path reports the same failure through its nested C calls; keep
  // this ownership regression focused on the context error and destruction.
  SUNContext_ClearErrHandlers(ctx);
  TestVector::fail_after = 2;
  N_Vector* clones = N_VCloneVectorArray(4, user.GetNVector());
  EXPECT_EQ(clones, nullptr);
  EXPECT_EQ(TestVector::destroyed, 2);
  EXPECT_EQ(TestVector::live, 1);
  EXPECT_EQ(SUNContext_GetLastError(ctx), SUN_ERR_OP_FAIL);
}

TEST_F(NVectorBaseTest, ExceptionFromErrOperationReachesHandler)
{
  TestVector user(ctx);
  sunrealtype c = 1;
  N_Vector input = user.GetNVector();
  N_Vector vectors[] = {input};
  EXPECT_EQ(input->ops->nvlinearcombination(1, &c, vectors, input), SUN_ERR_OP_FAIL);
  EXPECT_EQ(SUNContext_GetLastError(ctx), SUN_ERR_OP_FAIL);
  EXPECT_EQ(error.code, SUN_ERR_OP_FAIL);
  EXPECT_EQ(error.func, "N_VLinearCombination");
  EXPECT_EQ(error.message, "bridge test exception");
}
} // namespace
