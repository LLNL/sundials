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

#include <iostream>
#include <string>
#include <sundials/sundials_core.hpp>
#include <nvector/nvector_serial.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#define TTYPE int
#include "sundials_arraylist.h"

class SUNArrayListPODTest : public testing::Test {
 protected:
  SUNArrayList_int list;

  virtual void SetUp() override {
    list = SUNArrayList_int_New(2);
  }

  virtual void TearDown() override {
    SUNArrayList_int_Destroy(&list);
  }
};

TEST_F(SUNArrayListPODTest, NewAndDestroy) {
  EXPECT_NE(list, nullptr);
  EXPECT_EQ(list->size, 0);
  EXPECT_EQ(list->capacity, 2);
  SUNArrayList_int_Destroy(&list);
  EXPECT_EQ(list, nullptr);
}

TEST_F(SUNArrayListPODTest, IsEmpty) {
  EXPECT_TRUE(SUNArrayList_int_IsEmpty(list));
}

TEST_F(SUNArrayListPODTest, PushBack) {
  int value1 = 10;
  int value2 = 20;

  SUNArrayList_int_PushBack(list, value1);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNArrayList_int_At(list, 0), value1);

  SUNArrayList_int_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);
  EXPECT_EQ(*SUNArrayList_int_At(list, 1), value2);

  // Test resize
  int value3 = 30;
  SUNArrayList_int_PushBack(list, value3);
  EXPECT_EQ(list->size, 3);
  EXPECT_EQ(list->capacity, 3);
  EXPECT_EQ(*SUNArrayList_int_At(list, 2), value3);
}

TEST_F(SUNArrayListPODTest, At) {
  int value = 10;
  SUNArrayList_int_PushBack(list, value);

  EXPECT_EQ(*SUNArrayList_int_At(list, 0), value);
  // Test out of bounds
  EXPECT_EQ(SUNArrayList_int_At(list, -1), nullptr);
  EXPECT_EQ(SUNArrayList_int_At(list, 2), nullptr);
}

TEST_F(SUNArrayListPODTest, Set) {
  int value1 = 10;
  int value2 = 20;
  SUNArrayList_int_PushBack(list, value1);
  SUNArrayList_int_PushBack(list, value2);

  EXPECT_EQ(*SUNArrayList_int_At(list, 0), value1);
  SUNArrayList_int_Set(list, 0, value2);
  EXPECT_EQ(*SUNArrayList_int_At(list, 0), value2);
  // Test out of bounds
  SUNArrayList_int_Set(list, -1, value1);  // No effect
  SUNArrayList_int_Set(list, 2, value1);  // No effect
}

TEST_F(SUNArrayListPODTest, PopBack) {
  int value1 = 10;
  int value2 = 20;
  SUNArrayList_int_PushBack(list, value1);
  SUNArrayList_int_PushBack(list, value2);

  EXPECT_EQ(list->size, 2);
  SUNArrayList_int_PopBack(list);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNArrayList_int_At(list, 0), value1);
  // Pop from empty list
  SUNArrayList_int_PopBack(list);
  EXPECT_EQ(list->size, 0);
}

#define TTYPE N_Vector
#include "sundials_arraylist.h"

class SUNArrayListComplexTest : public testing::Test {
 protected:
  SUNArrayList_N_Vector list;
  sundials::Context sunctx;

  virtual void SetUp() override {
    list = SUNArrayList_N_Vector_New(2);
  }

  virtual void TearDown() override {
    SUNArrayList_N_Vector_Destroy(&list);
  }
};

TEST_F(SUNArrayListComplexTest, NewAndDestroy) {
  EXPECT_NE(list, nullptr);
  EXPECT_EQ(list->size, 0);
  EXPECT_EQ(list->capacity, 2);
  SUNArrayList_N_Vector_Destroy(&list);
  EXPECT_EQ(list, nullptr);
}

TEST_F(SUNArrayListComplexTest, IsEmpty) {
  EXPECT_TRUE(SUNArrayList_N_Vector_IsEmpty(list));
}

TEST_F(SUNArrayListComplexTest, PushBack) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNArrayList_N_Vector_PushBack(list, value1);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 0), value1);
  EXPECT_EQ(N_VGetLength(*SUNArrayList_N_Vector_At(list, 0)), N_VGetLength(value1));

  SUNArrayList_N_Vector_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 1), value2);
  EXPECT_EQ(N_VGetLength(*SUNArrayList_N_Vector_At(list, 1)), N_VGetLength(value2));

  // Test resize
  N_Vector value3 = N_VNew_Serial(3, sunctx);
  SUNArrayList_N_Vector_PushBack(list, value3);
  EXPECT_EQ(list->size, 3);
  EXPECT_EQ(list->capacity, 3);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 2), value3);
  EXPECT_EQ(N_VGetLength(*SUNArrayList_N_Vector_At(list, 2)), N_VGetLength(value3));

  N_VDestroy(value1);
  N_VDestroy(value2);
  N_VDestroy(value3);
}

TEST_F(SUNArrayListComplexTest, At) {
  N_Vector value = N_VNew_Serial(1, sunctx);

  SUNArrayList_N_Vector_PushBack(list, value);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 0), value);

  // Test out of bounds
  EXPECT_EQ(SUNArrayList_N_Vector_At(list, -1), nullptr);
  EXPECT_EQ(SUNArrayList_N_Vector_At(list, 2), nullptr);

  N_VDestroy(value);
}

TEST_F(SUNArrayListComplexTest, Set) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNArrayList_N_Vector_PushBack(list, value1);
  SUNArrayList_N_Vector_PushBack(list, value2);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 0), value1);
  SUNArrayList_N_Vector_Set(list, 0, value2);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 0), value2);

  // Test out of bounds
  SUNArrayList_N_Vector_Set(list, -1, value1);
  SUNArrayList_N_Vector_Set(list, 2, value1);

  N_VDestroy(value1);
  N_VDestroy(value2);
}

TEST_F(SUNArrayListComplexTest, PopBack) {
  N_Vector value1 = N_VNew_Serial(1, sunctx);
  N_Vector value2 = N_VNew_Serial(2, sunctx);

  SUNArrayList_N_Vector_PushBack(list, value1);
  SUNArrayList_N_Vector_PushBack(list, value2);
  EXPECT_EQ(list->size, 2);

  SUNArrayList_N_Vector_PopBack(list);
  EXPECT_EQ(list->size, 1);
  EXPECT_EQ(*SUNArrayList_N_Vector_At(list, 0), value1);

  // Pop from empty list
  SUNArrayList_N_Vector_PopBack(list);
  EXPECT_EQ(list->size, 0);

  N_VDestroy(value1);
  N_VDestroy(value2);
}
