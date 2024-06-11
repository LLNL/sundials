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
#include <string>
#include <sundials/sundials_core.h>
#include <sunadjoint/sunadjoint_checkpointscheme_basic.h>
#include <nvector/nvector_serial.h>

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
    SUNContext_Free(&sunctx);
    N_VDestroy(state);
  }

  SUNContext sunctx;
  N_Vector state;
};

TEST_F(SUNAdjointCheckpointSchemeBasic, CreateWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1,
                                                100, SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  // err = SUNAdjointCheckpointScheme_Destroy(&cs);
  // EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertSingleStageWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1,
                                                100, SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunrealtype* state_data = N_VGetArrayPointer(state);
  N_VConst(sunrealtype{1.0}, state);

  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, -1, 0, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 0, 0, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  N_Vector loaded_vec = NULL;
  SUNAdjointCheckpointScheme_LoadVector(cs, 0, 0, &loaded_vec);

  sunrealtype* loaded_data = N_VGetArrayPointer(loaded_vec);
  for (sunindextype i = 0; i < N_VGetLength(state); ++i) {
    EXPECT_TRUE(loaded_data[i] == state_data[i]);
  }

  // err = SUNAdjointCheckpointScheme_Destroy(&cs);
  // EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNAdjointCheckpointSchemeBasic, InsertTwoStagesWorks)
{
  SUNErrCode err;
  SUNAdjointCheckpointScheme cs = NULL;

  err = SUNAdjointCheckpointScheme_Create_Basic(SUNDATAIOMODE_INMEM, 1,
                                                100, SUNTRUE, SUNFALSE, sunctx, &cs);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunrealtype* state_data = N_VGetArrayPointer(state);
  N_VConst(sunrealtype{1.0}, state);

  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, -1, 0, state);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 0, 0, state);
  EXPECT_EQ(err, SUN_SUCCESS);
  err = SUNAdjointCheckpointScheme_InsertVector(cs, 0, 1, 0, state);
  EXPECT_EQ(err, SUN_SUCCESS);

  N_Vector loaded_vec = NULL;
  SUNAdjointCheckpointScheme_LoadVector(cs, 0, 1, &loaded_vec);

  sunrealtype* loaded_data = N_VGetArrayPointer(loaded_vec);
  for (sunindextype i = 0; i < N_VGetLength(state); ++i) {
    EXPECT_TRUE(loaded_data[i] == state_data[i]);
  }

  // err = SUNAdjointCheckpointScheme_Destroy(&cs);
  // EXPECT_EQ(err, SUN_SUCCESS);
}