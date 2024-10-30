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
#include <sys/types.h>

#include <nvector/nvector_serial.h>
#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sunadjoint/sunadjoint_checkpointscheme_fixed.h>
#include <sundials/sundials_core.h>
#include <sunmemory/sunmemory_system.h>

static bool compare_vectors(N_Vector expected, N_Vector actual)
{
  sunrealtype* adata = N_VGetArrayPointer(actual);
  sunrealtype* edata = N_VGetArrayPointer(expected);
  for (sunindextype i = 0; i < N_VGetLength(expected); ++i)
  {
    if (edata[i] != adata[i])
    {
      fprintf(stderr, "compare_vectors\nexpected:\n");
      N_VPrint(expected);
      fprintf(stderr, "compare_vectors\nactual:\n");
      N_VPrint(actual);
      return false;
    }
  }
  return true;
}

static void fake_mutlistage_method(SUNAdjointCheckpointScheme cs, int steps,
                                   int stages, bool test_load = false)
{
  N_Vector state  = N_VNew_Serial(10, cs->sunctx);
  N_Vector loaded = N_VClone(state);

  sunrealtype t    = 0.0;
  sunrealtype tout = 0.0;
  sunrealtype dt   = 0.1;

  // Initial condition
  SUNErrCode err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 0, t, state);

  // Fake a multistage method checkpointing pattern
  for (int step = 0; step < steps; ++step)
  {
    for (int stage = 1; stage <= stages; ++stage)
    {
      N_VConst(step * stage, state);
      sunrealtype ts = t; //t + dt / stages;
      err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, ts, state);
      EXPECT_EQ(err, SUN_SUCCESS);
    }

    int stage_idx = step == 0 ? stages + 1 : stages;
    N_VConst(step * stage_idx, state);
    err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage_idx, t + dt,
                                                  state);
    EXPECT_EQ(err, SUN_SUCCESS);

    t += dt;
  }

  if (test_load)
  {
    t = dt * steps;
    for (int step = steps - 1; step >= 0; --step)
    {
      int stage_idx = step == 0 ? stages + 1 : stages;
      N_VConst(step * stage_idx, state);
      err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage_idx, 0,
                                                  &loaded, &tout);
      EXPECT_EQ(err, SUN_SUCCESS);
      EXPECT_EQ(t, tout);
      EXPECT_TRUE(compare_vectors(state, loaded));

      for (int stage = stages; stage >= 1; --stage)
      {
        N_VConst(step * stage, state);
        stage_idx = step == 0 ? stage : stage - 1;
        err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage_idx, 0,
                                                    &loaded, &tout);
        EXPECT_EQ(err, SUN_SUCCESS);
        EXPECT_EQ(t - dt, tout);
        EXPECT_TRUE(compare_vectors(state, loaded));
      }

      t -= dt;
    }
  }

  N_VDestroy(state);
  N_VDestroy(loaded);
}

class SUNAdjointCheckpointSchemeBasic : public testing::Test
{
protected:
  SUNAdjointCheckpointSchemeBasic()
  {
    SUNContext_Create(SUN_COMM_NULL, &sunctx);
    state        = N_VNew_Serial(10, sunctx);
    loaded_state = N_VClone(state);
    mem_helper   = SUNMemoryHelper_Sys(sunctx);
  }

  ~SUNAdjointCheckpointSchemeBasic()
  {
    N_VDestroy(state);
    N_VDestroy(loaded_state);
    SUNMemoryHelper_Destroy(mem_helper);
    SUNContext_Free(&sunctx);
  }

  SUNContext sunctx;
  SUNMemoryHelper mem_helper;
  N_Vector state;
  N_Vector loaded_state;
};

TEST_F(SUNAdjointCheckpointSchemeBasic, CreateWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  int64_t interval                  = 1;
  int64_t estimate                  = 1;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, SingleStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 10;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 1, 1, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 1, 2, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStepsWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 2, 1, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStepsTwoStagesWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 2, 2, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, SingleStageWithDeleteWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 1, 1, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStagesWithDeleteWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 1, 2, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStepsWithDeleteWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 2, 1, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, TwoStepsTwoStagesWithDeleteWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 2, 2, true);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, CanStillInsertAfterDeleting)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  sunrealtype tout                  = 0.0;
  int64_t interval                  = 1;
  int64_t estimate                  = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  fake_mutlistage_method(cs, 2, 1, false);

  // Load the last step
  int64_t step  = 1;
  int64_t stage = 1;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, 0, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step again
  tout = 10.0;
  N_VConst(sunrealtype{10.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, tout, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Load it again
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, 0, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(10.0, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}
