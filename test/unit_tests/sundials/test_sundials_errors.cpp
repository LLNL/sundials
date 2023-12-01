/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

class SUNErrHandlerFnTest : public testing::Test
{
protected:
  SUNErrHandlerFnTest() { SUNContext_Create(SUN_COMM_NULL, &sunctx); }

  ~SUNErrHandlerFnTest() { SUNContext_Free(&sunctx); }

  SUNContext sunctx;
};

TEST_F(SUNErrHandlerFnTest, SUNLogErrHandlerFnLogsWhenCalled)
{
  testing::internal::CaptureStderr();
  std::string message = "Test log handler";
  SUNLogErrHandlerFn(__LINE__, __func__, __FILE__, message.c_str(), -1, nullptr,
                     sunctx);
  std::string output = testing::internal::GetCapturedStderr();
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(__func__),
                                     testing::HasSubstr("Test log handler")));
}

TEST_F(SUNErrHandlerFnTest, SUNAbortErrHandlerFnAbortsWhenCalled)
{
  ASSERT_DEATH(
    {
      SUNAbortErrHandlerFn(__LINE__, __func__, __FILE__, "Test abort handler",
                           -1, nullptr, sunctx);
    },
    "SUNAbortErrHandler: Calling abort now, use a different error handler to "
    "avoid program termination.\n");
}

TEST_F(SUNErrHandlerFnTest, SUNAssertErrHandlerFnAbortsWhenCalled)
{
  ASSERT_DEATH(
    {
      SUNAssertErrHandlerFn(__LINE__, __func__, __FILE__, "Test assert handler",
                            -1, nullptr, sunctx);
    },
    "SUNAssertErrHandler: assert(.*) failed... terminating\n");
}

class SUNContextErrFunctionTests : public testing::Test
{
protected:
  SUNContextErrFunctionTests() { SUNContext_Create(SUN_COMM_NULL, &sunctx); }

  ~SUNContextErrFunctionTests() { SUNContext_Free(&sunctx); }

  SUNContext sunctx;
};

void firstHandler(int line, const char* func, const char* file, const char* msg,
                  SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  std::vector<int>* order = static_cast<std::vector<int>*>(err_user_data);
  order->push_back(0);
}

void secondHandler(int line, const char* func, const char* file, const char* msg,
                   SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  std::vector<int>* order = static_cast<std::vector<int>*>(err_user_data);
  order->push_back(1);
}

void thirdHandler(int line, const char* func, const char* file, const char* msg,
                  SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  std::vector<int>* order = static_cast<std::vector<int>*>(err_user_data);
  order->push_back(2);
}

TEST_F(SUNContextErrFunctionTests, SUNContextPushErrHandlerWorks)
{
  std::vector<int> order = {};
  SUNContext_ClearErrHandlers(sunctx);
  SUNContext_PushErrHandler(sunctx, firstHandler, static_cast<void*>(&order));
  SUNContext_PushErrHandler(sunctx, secondHandler, static_cast<void*>(&order));
  SUNContext_PushErrHandler(sunctx, thirdHandler, static_cast<void*>(&order));
  SUNHandleErr(__LINE__, __func__, __FILE__, -1, sunctx);
  EXPECT_EQ(order.size(), 3);
  EXPECT_EQ(order.at(0), 2);
  EXPECT_EQ(order.at(1), 1);
  EXPECT_EQ(order.at(2), 0);
}

TEST_F(SUNContextErrFunctionTests, SUNContextPopErrHandlerWorks)
{
  std::vector<int> order = {};
  SUNContext_ClearErrHandlers(sunctx);
  SUNContext_PushErrHandler(sunctx, firstHandler, static_cast<void*>(&order));
  SUNContext_PushErrHandler(sunctx, secondHandler, static_cast<void*>(&order));
  SUNContext_PushErrHandler(sunctx, thirdHandler, static_cast<void*>(&order));
  SUNContext_PopErrHandler(sunctx);
  SUNHandleErr(__LINE__, __func__, __FILE__, -1, sunctx);
  EXPECT_EQ(order.size(), 2);
  EXPECT_EQ(order.at(0), 1);
  EXPECT_EQ(order.at(1), 0);
}
