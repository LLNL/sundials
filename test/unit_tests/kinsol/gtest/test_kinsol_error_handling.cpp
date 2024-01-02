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
#include <kinsol/kinsol.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "kinsol/kinsol_impl.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_error_handling.err"};

class KINErrConditionTest : public testing::Test
{
protected:
  KINErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    kinmem = KINCreate(sunctx);
  }

  ~KINErrConditionTest() { KINFree(&kinmem); }

  void* kinmem;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(KINErrConditionTest, WarningIsPrinted)
{
  SUNLogger_SetWarningFilename(logger, errfile.c_str());
  KINMemRec* kin_mem = (KINMemRec*)kinmem;
  KINProcessError(kin_mem, KIN_WARNING, __LINE__, __func__, __FILE__, "test");
  SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
}

TEST_F(KINErrConditionTest, ErrorIsPrinted)
{
  SUNLogger_SetErrorFilename(logger, errfile.c_str());
  // -1 is an illegal value
  KINSetNumMaxIters(kinmem, -1);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSG_BAD_MXITER)));
}
