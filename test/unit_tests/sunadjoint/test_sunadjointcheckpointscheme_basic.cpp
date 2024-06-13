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

#include "sundials/sundials_errors.h"

bool compare_vectors(N_Vector expected, N_Vector actual)
{
  sunrealtype* adata = N_VGetArrayPointer(actual);
  sunrealtype* edata = N_VGetArrayPointer(expected);
  for (sunindextype i = 0; i < N_VGetLength(expected); ++i)
  {
    if (edata[i] != adata[i]) return false;
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
  }

  ~SUNAdjointCheckpointSchemeBasic()
  {
    N_VDestroy(state);
    N_VDestroy(loaded_state);
    SUNContext_Free(&sunctx);
  }

  SUNContext sunctx;
  N_Vector state;
  N_Vector loaded_state;
};

TEST_F(SUNAdjointCheckpointSchemeBasic, CreateWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNTRUE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertSingleStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNTRUE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  sunindextype step  = 0;
  sunindextype stage = -1;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert first stage solution
  stage++;
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNTRUE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  sunindextype step  = 0;
  sunindextype stage = -1;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert first stage solution
  stage++;
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert second stage solution
  stage++;
  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStepsWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNTRUE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  sunindextype step  = 0;
  sunindextype stage = -1;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert first stage solution
  stage++;
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert the second step solution
  step++;
  stage = -1;
  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Insert first stage solution
  stage++;
  N_VConst(sunrealtype{3.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, AreDeletedWhenNotKeeping)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  sunindextype step  = 0;
  sunindextype stage = -1;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Try to load it again, it should be deleted
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_ERR_CHECKPOINT_NOT_FOUND);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, CanStillInsertAfterDeleting)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert the step solution
  sunindextype step  = 0;
  sunindextype stage = -1;
  N_VConst(sunrealtype{0.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Insert first stage solution
  stage++;
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  // Try to load it again, it should be deleted
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_ERR_CHECKPOINT_NOT_FOUND);

  // Insert the second step solution
  step++;
  stage = -1;
  N_VConst(sunrealtype{2.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, step, stage, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, step, stage, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  EXPECT_TRUE(compare_vectors(state, loaded_state));

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}
