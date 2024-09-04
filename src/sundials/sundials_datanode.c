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

#include "sundatanode/sundatanode_inmem.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_memory.h"
#include "sundials_datanode.h"

SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode self;
  self = (SUNDataNode)malloc(sizeof(*self));
  SUNAssert(self, SUN_ERR_MEM_FAIL);

  SUNDataNode_Ops ops;
  ops = (SUNDataNode_Ops)malloc(sizeof(*ops));
  SUNAssert(self, SUN_ERR_MEM_FAIL);

  ops->hasChildren = NULL;
  ops->isLeaf      = NULL;
  ops->isList      = NULL;
  ops->isObject    = NULL;
  ops->addChild    = NULL;
  ops->getChild    = NULL;
  ops->removeChild = NULL;
  ops->getData     = NULL;
  ops->setData     = NULL;
  ops->destroy     = NULL;

  self->dtype   = 0;
  self->ops     = ops;
  self->content = NULL;
  self->sunctx  = sunctx;

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

  if (self->ops->isLeaf) { return self->ops->isLeaf(self, yes_or_no); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_IsList(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->isList) { return self->ops->isList(self, yes_or_no); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_HasChildren(const SUNDataNode self,
                                   sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->hasChildren)
  {
    return self->ops->hasChildren(self, yes_or_no);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddChild(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->addChild) { return self->ops->addChild(self, child_node); }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode self, const char* name,
                                     SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->addNamedChild)
  {
    return self->ops->addNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetChild(const SUNDataNode self, sundataindex_t index,
                                SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->getChild)
  {
    return self->ops->getChild(self, index, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode self, const char* name,
                                     SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->getNamedChild)
  {
    return self->ops->getNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode self, const char* name,
                                        SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->removeNamedChild)
  {
    return self->ops->removeNamedChild(self, name, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveChild(SUNDataNode self, sundataindex_t index,
                                   SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->removeChild)
  {
    return self->ops->removeChild(self, index, child_node);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetData(const SUNDataNode self, void** data,
                               size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->getData)
  {
    return self->ops->getData(self, data, data_stride, data_bytes);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetDataNvector(SUNDataNode self, N_Vector v, sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->getDataNvector)
  {
    return self->ops->getDataNvector(self, v, t);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetData(SUNDataNode self, SUNMemoryType src_mem_type,
                               SUNMemoryType node_mem_type, void* data,
                               size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->setData)
  {
    return self->ops->setData(self, src_mem_type, node_mem_type, data,
                              data_stride, data_bytes);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetDataNvector(SUNDataNode self, N_Vector v, sunrealtype t)
{
  SUNFunctionBegin(self->sunctx);

  if (self->ops->setDataNvector)
  {
    return self->ops->setDataNvector(self, v, t);
  }

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  if ((*node)->ops->destroy) { return (*node)->ops->destroy(node); }

  free(*node);
  *node = NULL;

  return SUN_SUCCESS;
}
