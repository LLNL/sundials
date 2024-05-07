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

#include <arkode/arkode.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_nvector.h"

static const std::string errfile{"test_arkode_error_handling.err"};

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  return 0;
}

class ARKodeErrConditionTest : public testing::Test
{
protected:
  ARKodeErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    v          = N_VNew_Serial(1, sunctx);
    arkode_mem = ARKStepCreate(f, f, sunrealtype{0.0}, v, sunctx);
  }

  ~ARKodeErrConditionTest()
  {
    N_VDestroy(v);
    ARKodeFree(&arkode_mem);
  }

  void* arkode_mem;
  N_Vector v;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(ARKodeErrConditionTest, WarningIsPrinted)
{
  SUNLogger_SetWarningFilename(logger, errfile.c_str());
  ARKodeMemRec* ark_mem = (ARKodeMemRec*)arkode_mem;
  arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__, "test");
  SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
}

TEST_F(ARKodeErrConditionTest, ErrorIsPrinted)
{
  SUNLogger_SetErrorFilename(logger, errfile.c_str());
  // negative reltol is illegal
  ARKodeSStolerances(arkode_mem, /* reltol= */ -1e-4, /* abstol= */ 1e-4);
  std::string output = dumpstderr(sunctx, errfile);
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSG_ARK_BAD_RELTOL)));
}
