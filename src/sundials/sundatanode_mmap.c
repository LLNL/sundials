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

#include "sundatanode_mmap.h"

#define GET_IMPL(node) ((SUNDataNode_MmapImpl) (node)->impl)
#define GET_PROP(node, prop) (GET_IMPL(node)->prop)

static SUNDataNode sunDataNodeMmap_CreateEmpty(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCallNoRet(SUNDataNode_CreateEmpty(sunctx, &node));

  node->hasChildren = SUNDataNode_HasChildren_Mmap;
  node->isLeaf = SUNDataNode_IsLeaf_Mmap;
  node->isList = SUNDataNode_IsList_Mmap;
  node->isObject = NULL;
  node->addChild = SUNDataNode_AddChild_Mmap;
  node->getChild = SUNDataNode_GetChild_Mmap;
  node->removeChild = SUNDataNode_RemoveChild_Mmap;
  node->getData = SUNDataNode_GetData_Mmap;
  node->setData = SUNDataNode_SetData_Mmap;
  node->destroy = SUNDataNode_Destroy_Mmap;

  SUNDataNode_MmapImpl impl = (SUNDataNode_MmapImpl)malloc(sizeof(struct SUNDataNode_MmapImpl_s));
  SUNAssertNoRet(impl, SUN_ERR_MEM_FAIL);

  node->impl = (void*)impl;

  return node;
}

SUNErrCode SUNDataNode_CreateList_Mmap(sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeMmap_CreateEmpty(sunctx);

  GET_PROP(node, anon_children) = (SUNDataNode*)malloc(sizeof(*node) * num_elements);
  SUNAssert(GET_PROP(node, anon_children), SUN_ERR_MEM_FAIL);
  GET_PROP(node, max_anon_children) = num_elements;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf_Mmap(void* leaf_data, size_t data_stride, size_t data_bytes,
                                       SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeMmap_CreateEmpty(sunctx);

  GET_PROP(node, leaf_data) = leaf_data;
  GET_PROP(node, data_stride) = data_stride;
  GET_PROP(node, data_bytes) = data_bytes;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  *yes_or_no = SUNFALSE;
  if (GET_PROP(node, leaf_data)) {
    *yes_or_no = SUNTRUE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsList_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  *yes_or_no = SUNFALSE;
  if (GET_PROP(node, anon_children)) {
    *yes_or_no = SUNTRUE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_HasChildren_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);
  *yes_or_no = GET_PROP(node, num_anon_children) != 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddChild_Mmap(SUNDataNode node, SUNDataNode child_node)
{
  SUNFunctionBegin(node->sunctx);

  sunbooleantype is_leaf;
  SUNCheckCall(SUNDataNode_IsLeaf_Mmap(node, &is_leaf));

  if (is_leaf) {
    return SUN_ERR_DATANODE_NODEISLEAF;
  }

  if (GET_PROP(node, num_anon_children) == GET_PROP(node, max_anon_children)) {
    return SUN_ERR_DATANODE_MAXCHILDREN;
  }

  GET_PROP(node, anon_children)[GET_PROP(node, num_anon_children)++] = child_node;
  GET_PROP(child_node, parent) = node;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetChild_Mmap(const SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_Mmap(node, &has_children));

  if (has_children) {
    *child_node = GET_PROP(node, anon_children)[index];
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_RemoveChild_Mmap(SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_Mmap(node, &has_children));

  if (has_children) {
    *child_node = GET_PROP(node, anon_children)[index];
    GET_PROP(*child_node, parent) = NULL;
    GET_PROP(node, anon_children)[index] = NULL;
    GET_PROP(node, num_anon_children)--;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetData_Mmap(const SUNDataNode node, void** data)
{
  SUNFunctionBegin(node->sunctx);

  *data = GET_PROP(node, leaf_data);

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetData_Mmap(SUNDataNode node, void* data, size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(node->sunctx);

  sunbooleantype is_list;
  SUNCheckCall(SUNDataNode_IsList_Mmap(node, &is_list));

  if (is_list) {
    return SUN_ERR_DATANODE_NODEISLIST;
  }

  GET_PROP(node, leaf_data) = data;
  GET_PROP(node, data_stride) = data_stride;
  GET_PROP(node, data_bytes) = data_bytes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_Destroy_Mmap(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  free(GET_PROP(*node, anon_children));

  return SUN_SUCCESS;
}
