/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
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

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateLeaf_InMem(mem_helper, sunctx, node_out));
    break;
  default:
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode,
                                  sundataindex num_elements, SUNContext sunctx,
                                  SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateList_InMem(num_elements, sunctx, node_out));
    break;
  default:
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode,
                                    sundataindex num_elements,
                                    SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  switch (io_mode)
  {
  case (SUNDATAIOMODE_INMEM):
    SUNCheckCall(SUNDataNode_CreateObject_InMem(num_elements, sunctx, node_out));
    break;
  default:
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->isLeaf)
  {
    SUNErrCode err = self->ops->isLeaf(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);

  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_IsList(const SUNDataNode self, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->isList)
  {
    SUNErrCode err = self->ops->isList(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_HasChildren(const SUNDataNode self,
                                   sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->hasChildren)
  {
    SUNErrCode err = self->ops->hasChildren(self, yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddChild(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->addChild)
  {
    SUNErrCode err = self->ops->addChild(self, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode self, const char* name,
                                     SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->addNamedChild)
  {
    SUNErrCode err = self->ops->addNamedChild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetChild(const SUNDataNode self, sundataindex index,
                                SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getChild)
  {
    SUNErrCode err = self->ops->getChild(self, index, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode self, const char* name,
                                     SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getNamedChild)
  {
    SUNErrCode err = self->ops->getNamedChild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveNamedChild(const SUNDataNode self, const char* name,
                                        SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->removeNamedChild)
  {
    SUNErrCode err = self->ops->removeNamedChild(self, name, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_RemoveChild(SUNDataNode self, sundataindex index,
                                   SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->removeChild)
  {
    SUNErrCode err = self->ops->removeChild(self, index, child_node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetData(const SUNDataNode self, void** data,
                               size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getData)
  {
    SUNErrCode err = self->ops->getData(self, data, data_stride, data_bytes);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_GetDataNvector(SUNDataNode self, N_Vector v, sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->getDataNvector)
  {
    SUNErrCode err = self->ops->getDataNvector(self, v, t);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetData(SUNDataNode self, SUNMemoryType src_mem_type,
                               SUNMemoryType node_mem_type, void* data,
                               size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->setData)
  {
    SUNErrCode err = self->ops->setData(self, src_mem_type, node_mem_type, data,
                                        data_stride, data_bytes);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_SetDataNvector(SUNDataNode self, N_Vector v, sunrealtype t)
{
  SUNFunctionBegin(self->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->setDataNvector)
  {
    SUNErrCode err = self->ops->setDataNvector(self, v, t);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNDataNode_Destroy(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if ((*node)->ops->destroy)
  {
    SUNErrCode err = (*node)->ops->destroy(node);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }

  free((*node)->ops);
  free(*node);
  *node = NULL;

  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_SUCCESS;
}
