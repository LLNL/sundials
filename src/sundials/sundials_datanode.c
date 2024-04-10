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
#include <sundials/priv/sundials_errors_impl.h>

#include "sundials_datanode.h"
#include "sundatanode_mmap.h"

SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  node = (SUNDataNode)malloc(sizeof(*node));
  SUNAssert(node, SUN_ERR_MEM_FAIL);

  node->hasChildren = NULL;
  node->isLeaf = NULL;
  node->isList = NULL;
  node->isObject = NULL;
  node->addChild = NULL;
  node->getChild = NULL;
  node->removeChild = NULL;
  node->getData = NULL;
  node->setData = NULL;
  node->destroy = NULL;
  node->impl = NULL;
  node->sunctx = sunctx;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf(SUNDataIOMode io_mode, void* leaf_data, size_t data_stride, size_t data_bytes,
                                  SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch(io_mode)
  {
    case(SUNDATAIOMODE_MMAP):
      SUNCheckCall(SUNDataNode_CreateLeaf_Mmap(leaf_data, data_stride, data_bytes, sunctx, node_out));
      break;
    default:
      return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode, sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch(io_mode)
  {
    case(SUNDATAIOMODE_MMAP):
      SUNCheckCall(SUNDataNode_CreateList_Mmap(num_elements, sunctx, node_out));
      break;
    default:
      return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode, sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch(io_mode)
  {
    case(SUNDATAIOMODE_MMAP):
      SUNCheckCall(SUNDataNode_CreateObject_Mmap(num_elements, sunctx, node_out));
      break;
    default:
      return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  if (node->isLeaf) {
    return node->isLeaf(node, yes_or_no);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_IsList(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  if (node->isList) {
    return node->isList(node, yes_or_no);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_HasChildren(const SUNDataNode node, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(node->sunctx);

  if (node->hasChildren) {
    return node->hasChildren(node, yes_or_no);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddChild(SUNDataNode node, SUNDataNode child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->addChild) {
    return node->addChild(node, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode node, const char* name, SUNDataNode child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->addNamedChild) {
    return node->addNamedChild(node, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetChild(const SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->getChild) {
    return node->getChild(node, index, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode node, const char* name, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->getNamedChild) {
    return node->getNamedChild(node, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode node, const char* name, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->removeNamedChild) {
    return node->removeNamedChild(node, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveChild(SUNDataNode node, sundataindex_t index, SUNDataNode* child_node)
{
  SUNFunctionBegin(node->sunctx);

  if (node->removeChild) {
    return node->removeChild(node, index, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetData(const SUNDataNode node, void** data)
{
  SUNFunctionBegin(node->sunctx);

  if (node->getData) {
    return node->getData(node, data);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetData(SUNDataNode node, void* data, size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(node->sunctx);

  if (node->setData) {
    return node->setData(node, data, data_stride, data_bytes);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  if ((*node)->destroy) {
    SUNCheckCall((*node)->destroy(node));
  }

  free(*node);
  *node = NULL;

  return SUN_SUCCESS;
}
