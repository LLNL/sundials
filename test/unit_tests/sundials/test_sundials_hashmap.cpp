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

#include "sundials_hashmap_impl.h"

const size_t init_capacity = 2;

// Helper function to free memory for value
void FreeValue(void* ptr) {
  free(ptr);
}

class SUNHashMapTest : public testing::Test {
 protected:
  SUNHashMap map;

  virtual void SetUp() override {
    SUNHashMap_New(init_capacity, &map);
  }

  virtual void TearDown() override {
    SUNHashMap_Destroy(&map, FreeValue);
  }
};

TEST_F(SUNHashMapTest, CapacityWorks) {
  EXPECT_EQ(init_capacity, SUNHashMap_Capacity(map));
}

TEST_F(SUNHashMapTest, InsertAndGetWorks) {
  SUNErrCode err = SUN_SUCCESS;
  const char* key = "test_key";
  int value = 42;

  err = SUNHashMap_Insert(map, key, &value);
  ASSERT_EQ(err, SUN_SUCCESS);

  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key, &retrieved_value);
  ASSERT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(value, *((int*)retrieved_value));
}

TEST_F(SUNHashMapTest, InsertAndGetTwiceWorks) {
  SUNErrCode err = SUN_SUCCESS;
  const char* key1 = "test_key1";
  const char* key2 = "test_key2";
  int value1 = 42;
  int value2 = 43;

  err = SUNHashMap_Insert(map, key1, &value1);
  ASSERT_EQ(err, 0);

  err = SUNHashMap_Insert(map, key2, &value2);
  ASSERT_EQ(err, 0);

  // Ensure resize happened
  ASSERT_EQ(SUNHashMap_Capacity(map), 2);

  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key1, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value1, *((int*)retrieved_value));

  err = SUNHashMap_GetValue(map, key2, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value2, *((int*)retrieved_value));
}


TEST_F(SUNHashMapTest, InsertRequiringResizeWorks) {
  SUNErrCode err = SUN_SUCCESS;
  const char* key1 = "test_key1";
  const char* key2 = "test_key2";
  const char* key3 = "test_key3";
  int value1 = 42;
  int value2 = 43;
  int value3 = 44;

  err = SUNHashMap_Insert(map, key1, &value1);
  ASSERT_EQ(err, 0);

  err = SUNHashMap_Insert(map, key2, &value2);
  ASSERT_EQ(err, 0);

  // This should trigger a resize since init_capacity is 2
  err = SUNHashMap_Insert(map, key3, &value3);
  ASSERT_EQ(err, 0);

  // Ensure resize happened
  ASSERT_EQ(SUNHashMap_Capacity(map), 3);

  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key1, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value1, *((int*)retrieved_value));

  err = SUNHashMap_GetValue(map, key2, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value2, *((int*)retrieved_value));

  err = SUNHashMap_GetValue(map, key3, &retrieved_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(value3, *((int*)retrieved_value));
}

TEST_F(SUNHashMapTest, InsertDuplicateKeyFails) {
  int err;

  // Insert same key twice (should overwrite)
  const char* key = "test_key";
  int value1 = 42;
  int value2 = 100;

  err = SUNHashMap_Insert(map, key, &value1);
  ASSERT_EQ(err, 0);
  err = SUNHashMap_Insert(map, key, &value2);
  ASSERT_EQ(err, -2);
}

TEST_F(SUNHashMapTest, RemoveWorks) {
  int err;

  // Insert a key-value pair
  const char* key = "test_key";
  int value = 42;
  err = SUNHashMap_Insert(map, key, &value);
  ASSERT_EQ(err, 0);

  // Remove the key
  void* removed_value;
  err = SUNHashMap_Remove(map, key, &removed_value);
  ASSERT_EQ(err, 0);
  EXPECT_EQ(&value, removed_value);

  // Check if key is gone
  void* retrieved_value;
  err = SUNHashMap_GetValue(map, key, &retrieved_value);
  ASSERT_EQ(err, -2);
}

// Add more tests for other functionalities like iteration, sorting, etc.
