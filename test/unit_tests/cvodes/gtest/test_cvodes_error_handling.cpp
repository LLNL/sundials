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

#include "gmock/gmock.h"
#include <cvodes/cvodes.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "cvodes/cvodes_impl.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_error_handling.err"};

class CVodeErrConditionTest : public testing::Test
{
protected:
  CVodeErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
  }

  ~CVodeErrConditionTest() { CVodeFree(&cvode_mem); }

  void* cvode_mem;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(CVodeErrConditionTest, WarningIsPrinted)
{
  SUNLogger_SetWarningFilename(logger, errfile.c_str());
  CVodeMemRec* ark_mem = (CVodeMemRec*)cvode_mem;
  cvProcessError(ark_mem, CV_WARNING, __LINE__, __func__, __FILE__, "test");
  SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
}

TEST_F(CVodeErrConditionTest, ErrorIsPrinted)
{
  SUNLogger_SetErrorFilename(logger, errfile.c_str());
  // attempting to call CVodeSStolerances before CVodeInit is illegal
  CVodeSStolerances(cvode_mem, 1e-4, 1e-4);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSGCV_NO_MALLOC)));
}
