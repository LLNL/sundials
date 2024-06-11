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

#include <fstream>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <string>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <sundials/sundials_core.h>

void compare_vectors(N_Vector expected, N_Vector actual)
{
  sunrealtype* adata = N_VGetArrayPointer(actual);
  sunrealtype* edata = N_VGetArrayPointer(expected);
  for (sunindextype i = 0; i < N_VGetLength(expected); ++i)
  {
    EXPECT_TRUE(edata[i] == adata[i]);
  }
}

class SUNAdjointCheckpointSchemeBasic : public testing::Test
{
protected:
  SUNAdjointCheckpointSchemeBasic()
  {
    SUNContext_Create(SUN_COMM_NULL, &sunctx);
    state = N_VNew_Serial(10, sunctx);
  }

  ~SUNAdjointCheckpointSchemeBasic()
  {
    N_VDestroy(state);
    SUNContext_Free(&sunctx);
  }

  SUNContext sunctx;
  N_Vector state;
};

TEST_F(SUNAdjointCheckpointSchemeBasic, CreateWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNFALSE, sunctx, &cs);
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
                                                SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  N_VConst(sunrealtype{0.0}, state);
  N_Vector loaded_state = NULL;

  // Insert the step solution
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, -1, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, 0, -1, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  compare_vectors(state, loaded_state);

  // Insert first stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 0, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, 0, 0, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  compare_vectors(state, loaded_state);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStepsWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;
  const sunrealtype t           = 0.0;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1, 100,
                                                SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  N_VConst(sunrealtype{0.0}, state);
  N_Vector loaded_state = NULL;

  // Insert the step solution
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, -1, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, 0, -1, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  compare_vectors(state, loaded_state);

  // Insert first stage solution
  N_VConst(sunrealtype{1.0}, state);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 0, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, 0, 0, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  compare_vectors(state, loaded_state);

  // Insert the step solution
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 1, -1, t, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  // Try to load it
  err = SUNAdjointCheckpointScheme_LoadVector(cs, 0, 0, &loaded_state);
  EXPECT_EQ(err, SUN_SUCCESS);
  compare_vectors(state, loaded_state);

  err = SUNAdjointCheckpointScheme_Destroy(&cs);
  EXPECT_EQ(err, SUN_SUCCESS);
}