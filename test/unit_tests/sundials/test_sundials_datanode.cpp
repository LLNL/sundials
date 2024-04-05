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

  err = SUNDataNode_Destroy(&node);
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

  err = SUNDataNode_Destroy(&node);
  EXPECT_EQ(err, SUN_SUCCESS);
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

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
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

  EXPECT_EQ(root_node, child_node->parent);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, AddChildFailsWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 1;

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  child_node = NULL;
  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_ERR_DATANODE_NODEISLEAF);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, AddChildFailsCorrectlyWhenFull)
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

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}


TEST_F(SUNDataNodeTest, HasChildrenWorks)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunbooleantype yes_or_no = SUNTRUE;
  err = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_FALSE(yes_or_no);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_AddChild(root_node, child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_TRUE(yes_or_no);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}


TEST_F(SUNDataNodeTest, RemoveChildWorks)
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

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  sunbooleantype yes_or_no = SUNTRUE;
  err = SUNDataNode_HasChildren(root_node, &yes_or_no);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_FALSE(yes_or_no);

  EXPECT_FALSE(child_node->parent);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, RemoveChildWorksWhenEmpty)
{
  SUNErrCode err;
  SUNDataNode root_node, child_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_RemoveChild(root_node, 0, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetChildWorks)
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

  err = SUNDataNode_GetChild(root_node, 0, &child_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(integer_value, get_leaf_as_int(child_node));

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  err = SUNDataNode_Destroy(&child_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetDataWorksWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node;
  unsigned int num_elem = 1;

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  err = SUNDataNode_GetData(root_node, &raw_value);
  EXPECT_EQ(err, SUN_SUCCESS);

  int value_we_got = *((int*) raw_value);
  EXPECT_EQ(integer_value, value_we_got);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, GetDataWorksWhenList)
{
  SUNErrCode err;
  SUNDataNode root_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  err = SUNDataNode_GetData(root_node, &raw_value);
  EXPECT_EQ(err, SUN_SUCCESS);

  EXPECT_EQ(raw_value, (void*)NULL);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, SetDataFailsCorrectlyWhenList)
{
  SUNErrCode err;
  SUNDataNode root_node;
  unsigned int num_elem = 1;

  err = SUNDataNode_CreateList(num_elem, sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  err = SUNDataNode_SetData(root_node, &raw_value, 0, 0);
  EXPECT_EQ(err, SUN_ERR_DATANODE_NODEISLIST);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}

TEST_F(SUNDataNodeTest, SetDataWorksWhenLeaf)
{
  SUNErrCode err;
  SUNDataNode root_node;
  unsigned int num_elem = 1;

  int integer_value = 5;
  err = SUNDataNode_CreateLeaf((void*)&integer_value, 0, sizeof(integer_value), sunctx, &root_node);
  EXPECT_EQ(err, SUN_SUCCESS);

  int new_integer_value = 3;
  err = SUNDataNode_SetData(root_node, (void*)&new_integer_value, 0, sizeof(new_integer_value));
  EXPECT_EQ(err, SUN_SUCCESS);

  void* raw_value;
  err = SUNDataNode_GetData(root_node, &raw_value);
  EXPECT_EQ(err, SUN_SUCCESS);

  int value_we_got = *((int*) raw_value);
  EXPECT_EQ(new_integer_value, value_we_got);

  err = SUNDataNode_Destroy(&root_node);
  EXPECT_EQ(err, SUN_SUCCESS);
}
