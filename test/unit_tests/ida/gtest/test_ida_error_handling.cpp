/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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
#include <ida/ida.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "ida/ida_impl.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_error_handling.err"};

class IDAErrConditionTest : public testing::Test
{
protected:
  IDAErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    IDA_mem = IDACreate(sunctx);
  }

  ~IDAErrConditionTest() { IDAFree(&IDA_mem); }

  void* IDA_mem;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(IDAErrConditionTest, WarningIsPrinted)
{
  SUNLogger_SetWarningFilename(logger, errfile.c_str());
  IDAMemRec* ark_mem = (IDAMemRec*)IDA_mem;
  IDAProcessError(ark_mem, IDA_WARNING, __LINE__, __func__, __FILE__, "test");
  SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
}

TEST_F(IDAErrConditionTest, ErrorIsPrinted)
{
  SUNLogger_SetErrorFilename(logger, errfile.c_str());
  // attempting to call IDASStolerances before IDAInit is illegal
  IDASStolerances(IDA_mem, 1e-4, 1e-4);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSG_NO_MALLOC)));
}
