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
#include <sundials/sundials_datanode.h>
#include <nvector/nvector_serial.h>

#include "../utilities/dumpstderr.hpp"

int get_leaf_as_int(SUNDataNode node)
{
  return *((int*) node->leaf_data);
}

class SUNDataNodeTest : public testing::Test
{
protected:
  SUNDataNodeTest()
  {
    SUNContext_Create(SUN_COMM_NULL, &sunctx);
  }

  ~SUNDataNodeTest()
  {
    SUNContext_Free(&sunctx);
  }

  SUNContext sunctx;
};

TEST_F(SUNDataNodeTest, CreateEmptyWorks)
{
  SUNErrCode err;
  SUNDataNode node;

  err = SUNDataNode_CreateEmpty(sunctx, &node);
  EXPECT_EQ(err, SUN_SUCCESS);
}


TEST_F(SUNDataNodeTest, CreateLeafWorks)
{
  SUNErrCode err;
  SUNDataNode node;
  int integer_value = 5;

  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(integer_value, get_leaf_as_int(node));
}

TEST_F(SUNDataNodeTest, CreateListWorks)
{
  SUNErrCode err;
  SUNDataNode root_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(num_elem, root_node->max_anon_children);
  EXPECT_EQ(0, root_node->num_anon_children);
}

TEST_F(SUNDataNodeTest, AddChildWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(1, root_node->num_anon_children);

  EXPECT_EQ(integer_value, get_leaf_as_int(root_node->anon_children[0]));
}

TEST_F(SUNDataNodeTest, AddChildFailsCorrectly)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_ERR_DATANODE_MAXCHILDREN);

  EXPECT_EQ(1, root_node->num_anon_children);
  EXPECT_EQ(integer_value, get_leaf_as_int(root_node->anon_children[0]));
}

