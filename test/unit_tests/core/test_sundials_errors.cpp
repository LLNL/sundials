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

#include <gtest/gtest.h>
#include <sundials/sundials.h>

class SUNErrHandlerTest : public ::testing::Test
{
protected:
  SUNErrHandlerTest() { SUNContext_Create(nullptr, &sunctx); }

  ~SUNErrHandlerTest() { SUNContext_Free(&sunctx); }

  SUNContext sunctx;
};

TEST_F(SUNErrHandlerTest, SUNLogErrHandlerFnLogsWhenCalled)
{
  testing::internal::CaptureStderr();
  std::string expected = "[ERROR][rank 0][/Volumes/Workspace/SUNDIALS/repos/feature/error-handling/test/unit_tests/core/test_sundials_errors.cpp:30][TestBody] Test log handler\n";
  SUNLogErrHandlerFn(__LINE__, __func__, __FILE__, "Test log handler", -1,
                     nullptr, sunctx);
  std::string output = testing::internal::GetCapturedStderr(); 
  ASSERT_EQ(output, expected);
}

TEST_F(SUNErrHandlerTest, SUNAbortErrHandlerFnAbortsWhenCalled)
{
  ASSERT_DEATH(
    {
      SUNAbortErrHandlerFn(__LINE__, __func__, __FILE__, "Test abort handler",
                           -1, nullptr, sunctx);
    },
    "SUNAbortErrHandler: Calling abort now, use a different error handler to "
    "avoid program termination.\n");
}

TEST_F(SUNErrHandlerTest, SUNAssertErrHandlerFnAbortsWhenCalled)
{
  ASSERT_DEATH(
    {
      SUNAssertErrHandlerFn(__LINE__, __func__, __FILE__, "Test assert handler",
                            -1, nullptr, sunctx);
    },
    "SUNAssertErrHandler: assert(.*) failed... terminating\n");
}
