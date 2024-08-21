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
#include <nvector/nvector_serial.h>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sundials/sundials_core.h>
#include <sys/types.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_memory.h"
#include "sundials/sundials_types.h"
#include "sunmemory/sunmemory_system.h"

bool compare_vectors(N_Vector expected, N_Vector actual)
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

  uint64_t interval                 = 1;
  uint64_t estimate                 = 1;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertSingleStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  const sunrealtype t               = 1.0;
  uint64_t interval                 = 1;
  uint64_t estimate                 = 10;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the initial condition
  int64_t step  = 0;
  int64_t stage = 0;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  sunrealtype tout = 0.0;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  step++;

  // Insert first stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  const sunrealtype t               = 1.0;
  uint64_t interval                 = 1;
  uint64_t estimate                 = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the initial condition
  int64_t step  = 0;
  int64_t stage = 0;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  sunrealtype tout = 0.0;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  step++;

  // Insert first stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  stage++;

  // Insert second stage solution
  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(t, tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStepsWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  const sunrealtype t               = 0.0;
  uint64_t interval                 = 1;
  uint64_t estimate                 = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNTRUE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the initial condition
  int64_t step  = 0;
  int64_t stage = 0;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  sunrealtype tout = 0.0;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  step++;

  // Insert first stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  stage++;

  // Insert step solution
  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  step++;
  stage = 0;

  // Insert second step, first stage
  N_VConst(sunrealtype{3.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  stage++;

  // Insert second step solution
  N_VConst(sunrealtype{4.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, AreDeletedWhenNotKeeping)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  const sunrealtype t               = 0.0;
  uint64_t interval                 = 1;
  uint64_t estimate                 = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  int64_t step  = 0;
  int64_t stage = 0;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  sunrealtype tout = 0.0;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Try to load it again, it should be deleted
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(t, tout);
  EXPECT_EQ(err, SUN_ERR_CHECKPOINT_NOT_FOUND);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, CanStillInsertAfterDeleting)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs     = NULL;
  const sunrealtype t               = 0.0;
  uint64_t interval                 = 1;
  uint64_t estimate                 = 100;
  sunbooleantype save_stages        = SUNTRUE;
  sunbooleantype keep_after_loading = SUNFALSE;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, mem_helper,
                                                interval, estimate, save_stages,
                                                keep_after_loading, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the initial solution
  int64_t step  = 0;
  int64_t stage = 0;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  step++;

  // Insert first step, stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  sunrealtype tout = 0.0;
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Try to load it again, it should be deleted
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(t, tout);
  EXPECT_EQ(err, SUN_ERR_CHECKPOINT_NOT_FOUND);

  // Insert the second step solution
  step++;
  stage = 0;

  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state,
                                              &tout);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_EQ(t, tout);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}
