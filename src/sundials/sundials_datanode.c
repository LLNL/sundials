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

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_datanode.h>

#include "sundatanode/sundatanode_inmem.h"
#include "sundials/sundials_memory.h"

SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode self;
  self = (SUNDataNode)malloc(sizeof(*self));
  SUNAssert(self, SUN_ERR_MEM_FAIL);

  self->hasChildren = NULL;
  self->isLeaf      = NULL;
  self->isList      = NULL;
  self->isObject    = NULL;
  self->addChild    = NULL;
  self->getChild    = NULL;
  self->removeChild = NULL;
  self->getData     = NULL;
  self->setData     = NULL;
  self->destroy     = NULL;
  self->impl        = NULL;
  self->sunctx      = sunctx;

  *node_out = self;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf(SUNDataIOMode io_mode,
                                  SUNMemoryHelper mem_helper, SUNContext sunctx,
                                  SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateLeaf_InMem(mem_helper, sunctx, node_out));
    break;
  default: return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode,
                                  sundataindex_t num_elements,
                                  SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateList_InMem(num_elements, sunctx, node_out));
    break;
  default: return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode,
                                    sundataindex_t num_elements,
                                    SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateObject_InMem(num_elements, sunctx, node_out));
    break;
  default: return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (self->isLeaf) { return self->isLeaf(self, yes_or_no); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_IsList(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (self->isList) { return self->isList(self, yes_or_no); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_HasChildren(const SUNDataNode self,
                                   sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (self->hasChildren) { return self->hasChildren(self, yes_or_no); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddChild(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->addChild) { return self->addChild(self, child_node); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode self, const char* name,
                                     SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->addNamedChild)
  {
    return self->addNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetChild(const SUNDataNode self, sundataindex_t index,
                                SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->getChild) { return self->getChild(self, index, child_node); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode self, const char* name,
                                     SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->getNamedChild)
  {
    return self->getNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode self, const char* name,
                                        SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->removeNamedChild)
  {
    return self->removeNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveChild(SUNDataNode self, sundataindex_t index,
                                   SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->removeChild) { return self->removeChild(self, index, child_node); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetData(const SUNDataNode self, void** data,
                               size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  if (self->getData)
  {
    return self->getData(self, data, data_stride, data_bytes);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetDataNvector(SUNDataNode self, N_Vector v)
{
  SUNFunctionBegin(self->sunctx);

  if (self->getDataNvector) { return self->getDataNvector(self, v); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetData(SUNDataNode self, SUNMemoryType src_mem_type,
                               SUNMemoryType node_mem_type, void* data,
                               size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  if (self->setData)
  {
    return self->setData(self, src_mem_type, node_mem_type, data, data_stride,
                         data_bytes);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetDataNvector(SUNDataNode self, N_Vector v)
{
  SUNFunctionBegin(self->sunctx);

  if (self->setDataNvector) { return self->setDataNvector(self, v); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  if ((*node)->destroy) { return (*node)->destroy(node); }

  free(*node);
  *node = NULL;

  return SUN_SUCCESS;
}
