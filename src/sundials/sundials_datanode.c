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

#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>
#include <sundials/priv/sundials_errors_impl.h>


SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  node = (SUNDataNode)malloc(sizeof(*node));
  SUNAssert(node, SUN_ERR_MEM_FAIL);

  node->hasChildren = NULL;
  node->addChild = NULL;
  node->getChild = NULL;
  node->removeChild = NULL;
  node->getData = NULL;
  node->setData = NULL;
  node->destroy = NULL;

  node->parent = NULL;

  node->leaf_data = NULL;
  node->data_stride = 0;
  node->data_bytes = 0;

  // node->name = NULL;
  // node->named_children = NULL;
  // node->num_named_children = NULL;

  node->anon_children = NULL;
  node->max_anon_children = 0;
  node->num_anon_children = 0;

  node->impl = NULL;

  node->sunctx = sunctx;

  *node_out = node;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateList(sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCall(SUNDataNode_CreateEmpty(sunctx, &node));

  node->anon_children = (SUNDataNode*)malloc(sizeof(*node) * num_elements);
  SUNAssert(node, SUN_ERR_MEM_FAIL);

  node->max_anon_children = num_elements;

  *node_out = node;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf(void* leaf_data, size_t data_stride, size_t data_bytes, SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCall(SUNDataNode_CreateEmpty(sunctx, &node));

  node->leaf_data = leaf_data;
  node->data_stride = data_stride;
  node->data_bytes = data_bytes;

  *node_out = node;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_NodeIsLeaf(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  *yes_or_no = SUNFALSE;
  if (node->leaf_data) {
    *yes_or_no = SUNTRUE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_NodeIsList(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  *yes_or_no = SUNFALSE;
  if (node->anon_children) {
    *yes_or_no = SUNTRUE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_HasChildren(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);
  *yes_or_no = node->num_anon_children != 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddChild(SUNDataNode parent_node, SUNDataNode child_node)
{
  SUNFunctionBegin(parent_node->sunctx);

  sunbooleantype is_leaf;
  SUNCheckCall(SUNDataNode_NodeIsLeaf(parent_node, &is_leaf));

  if (is_leaf) {
    return SUN_ERR_DATANODE_NODEISLEAF;
  }

  if (parent_node->num_anon_children == parent_node->max_anon_children) {
    return SUN_ERR_DATANODE_MAXCHILDREN;
  }

  parent_node->anon_children[parent_node->num_anon_children++] = child_node;
  child_node->parent = parent_node;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetChild(const SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren(node, &has_children));

  if (has_children) {
    *child_node = node->anon_children[index];
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_RemoveChild(SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren(node, &has_children));

  if (has_children) {
    *child_node = node->anon_children[index];
    (*child_node)->parent = NULL;
    node->anon_children[index] = NULL;
    node->num_anon_children--;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetData(const SUNDataNode node, void** data)
{
  SUNFunctionBegin(node->sunctx);

  *data = node->leaf_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetData(SUNDataNode node, void* data, size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(node->sunctx);

  sunbooleantype is_list;
  SUNCheckCall(SUNDataNode_NodeIsList(node, &is_list));

  if (is_list) {
    return SUN_ERR_DATANODE_NODEISLIST;
  }

  node->leaf_data = data;
  node->data_stride = data_stride;
  node->data_bytes = data_bytes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  free((*node)->anon_children);
  free(*node);
  *node = NULL;

  return SUN_SUCCESS;
}
