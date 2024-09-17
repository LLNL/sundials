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

#include <string.h>

#include "sundatanode/sundatanode_inmem.h"
#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_memory.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sundials_datanode.h"
#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"

#define GET_IMPL(node)        ((SUNDataNode_InMemContent)(node)->content)
#define IMPL_PROP(node, prop) (GET_IMPL(node)->prop)
#define BASE_PROP(node, prop) ((node)->prop)

static SUNDataNode sunDataNodeInMem_CreateEmpty(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCallNoRet(SUNDataNode_CreateEmpty(sunctx, &node));

  node->ops->hasChildren      = SUNDataNode_HasChildren_InMem;
  node->ops->isLeaf           = SUNDataNode_IsLeaf_InMem;
  node->ops->isList           = SUNDataNode_IsList_InMem;
  node->ops->isObject         = SUNDataNode_IsObject_InMem;
  node->ops->addChild         = SUNDataNode_AddChild_InMem;
  node->ops->addNamedChild    = SUNDataNode_AddNamedChild_InMem;
  node->ops->getChild         = SUNDataNode_GetChild_InMem;
  node->ops->getNamedChild    = SUNDataNode_GetNamedChild_InMem;
  node->ops->removeChild      = SUNDataNode_RemoveChild_InMem;
  node->ops->removeNamedChild = SUNDataNode_RemoveNamedChild_InMem;
  node->ops->getData          = SUNDataNode_GetData_InMem;
  node->ops->getDataNvector   = SUNDataNode_GetDataNvector_InMem;
  node->ops->setData          = SUNDataNode_SetData_InMem;
  node->ops->setDataNvector   = SUNDataNode_SetDataNvector_InMem;
  node->ops->destroy          = SUNDataNode_Destroy_InMem;

  SUNDataNode_InMemContent content =
    (SUNDataNode_InMemContent)malloc(sizeof(struct SUNDataNode_InMemImpl_));
  SUNAssertNoRet(content, SUN_ERR_MEM_FAIL);

  content->parent             = NULL;
  content->mem_helper         = NULL;
  content->leaf_data          = NULL;
  content->name               = NULL;
  content->named_children     = NULL;
  content->num_named_children = 0;
  content->anon_children      = NULL;

  node->content = (void*)content;

  return node;
}

static void sunDataNodeInMem_DestroyEmpty(SUNDataNode* node)
{
  if (!node || !(*node)) { return; }
  if (BASE_PROP(*node, content)) { free(BASE_PROP(*node, content)); }
  BASE_PROP(*node, content) = NULL;
  free(*node);
  *node = NULL;
}

static void sunDataNodeFreeKeyValue(SUNHashMapKeyValue* kv_ptr);
static void sunDataNodeFreeValue(SUNDataNode* nodeptr);

SUNErrCode SUNDataNode_CreateList_InMem(sundataindex_t init_size,
                                        SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeInMem_CreateEmpty(sunctx);

  BASE_PROP(node, dtype) = SUNDATANODE_LIST;
  IMPL_PROP(node, anon_children) =
    SUNStlVector_SUNDataNode_New(init_size, sunDataNodeFreeValue);

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject_InMem(sundataindex_t init_size,
                                          SUNContext sunctx,
                                          SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeInMem_CreateEmpty(sunctx);

  BASE_PROP(node, dtype) = SUNDATANODE_OBJECT;

  SUNHashMap map;
  SUNCheckCall(SUNHashMap_New(init_size, sunDataNodeFreeKeyValue, &map));

  IMPL_PROP(node, named_children) = map;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf_InMem(SUNMemoryHelper mem_helper,
                                        SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeInMem_CreateEmpty(sunctx);

  BASE_PROP(node, dtype)      = SUNDATANODE_LEAF;
  IMPL_PROP(node, mem_helper) = mem_helper;
  IMPL_PROP(node, leaf_data)  = NULL;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsLeaf_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_PROP(self, dtype) == SUNDATANODE_LEAF;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsList_InMem(const SUNDataNode self,
                                    sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_PROP(self, dtype) == SUNDATANODE_LIST;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_IsObject_InMem(const SUNDataNode self,
                                      sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no = BASE_PROP(self, dtype) == SUNDATANODE_OBJECT;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_HasChildren_InMem(const SUNDataNode self,
                                         sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  *yes_or_no =
    (IMPL_PROP(self, anon_children) &&
     SUNStlVector_SUNDataNode_Size(IMPL_PROP(self, anon_children)) != 0) ||
    IMPL_PROP(self, num_named_children) != 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddChild_InMem(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_PROP(self, dtype) == SUNDATANODE_LIST, SUN_ERR_ARG_WRONGTYPE);
  SUNStlVector_SUNDataNode_PushBack(IMPL_PROP(self, anon_children), child_node);
  IMPL_PROP(child_node, parent) = self;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddNamedChild_InMem(SUNDataNode self, const char* name,
                                           SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_PROP(self, dtype) == SUNDATANODE_OBJECT, SUN_ERR_ARG_WRONGTYPE);

  IMPL_PROP(child_node, name) = name;
  if (SUNHashMap_Insert(IMPL_PROP(self, named_children), name, child_node))
  {
    // fprintf(stdout, "node with name=%s could not be inserted, current named children:\n",
    //         name);
    // SUNHashMap_PrintKeys(IMPL_PROP(self, named_children), stdout);
    return SUN_ERR_OP_FAIL;
  }

  IMPL_PROP(child_node, parent) = self;
  IMPL_PROP(self, num_named_children)++;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetChild_InMem(const SUNDataNode self,
                                      sundataindex_t index,
                                      SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (!has_children) { return SUN_ERR_DATANODE_NODENOTFOUND; }

  SUNDataNode* child_node_ptr =
    SUNStlVector_SUNDataNode_At(IMPL_PROP(self, anon_children), index);
  if (child_node_ptr) { *child_node = *child_node_ptr; }
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetNamedChild_InMem(const SUNDataNode self,
                                           const char* name,
                                           SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (has_children)
  {
    if (SUNHashMap_GetValue(IMPL_PROP(self, named_children), name,
                            (void**)child_node))
    {
      // fprintf(stdout, "node with name=%s not found, current named children:\n",
      //         name);
      // SUNHashMap_PrintKeys(IMPL_PROP(self, named_children), stdout);
      return SUN_ERR_DATANODE_NODENOTFOUND;
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_RemoveChild_InMem(SUNDataNode self, sundataindex_t index,
                                         SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (!has_children) { return SUN_SUCCESS; }

  SUNDataNode* child_node_ptr =
    SUNStlVector_SUNDataNode_At(IMPL_PROP(self, anon_children), index);
  if (child_node_ptr)
  {
    *child_node = *child_node_ptr;
    if (*child_node)
    {
      IMPL_PROP(*child_node, parent) = NULL;
      SUNStlVector_SUNDataNode_Erase(IMPL_PROP(self, anon_children), index);
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_RemoveNamedChild_InMem(const SUNDataNode self,
                                              const char* name,
                                              SUNDataNode* child_node)
{
  SUNFunctionBegin(self->sunctx);

  *child_node = NULL;

  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  if (has_children)
  {
    if (SUNHashMap_Remove(IMPL_PROP(self, named_children), name,
                          (void**)child_node))
    {
      return SUN_ERR_DATANODE_NODENOTFOUND;
    }
    IMPL_PROP(*child_node, parent) = NULL;
    IMPL_PROP(self, num_named_children)--;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetData_InMem(const SUNDataNode self, void** data,
                                     size_t* data_stride, size_t* data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNMemory leaf_data = (SUNMemory)IMPL_PROP(self, leaf_data);

  *data_stride = leaf_data->stride;
  *data_bytes  = leaf_data->bytes;
  *data        = leaf_data->ptr;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_GetDataNvector_InMem(const SUNDataNode self, N_Vector v,
                                            sunrealtype* t)
{
  SUNFunctionBegin(self->sunctx);

  void* queue = NULL;

  SUNMemory leaf_data = (SUNMemory)IMPL_PROP(self, leaf_data);

  SUNMemoryType leaf_mem_type = leaf_data->type;
  SUNMemoryType buffer_mem_type = N_VGetDeviceArrayPointer(v) ? SUNMEMTYPE_DEVICE
                                                              : SUNMEMTYPE_HOST;

  sunindextype buffer_size = 0;
  SUNCheckCall(N_VBufSize(v, &buffer_size));
  SUNAssert((buffer_size + sizeof(sunrealtype)) == leaf_data->bytes,
            SUN_ERR_ARG_INCOMPATIBLE);

  if (leaf_mem_type == buffer_mem_type)
  {
    sunrealtype* data_ptr = leaf_data->ptr;
    *t                    = data_ptr[0];
    SUNCheckCall(N_VBufUnpack(v, &data_ptr[1]));
  }
  else
  {
    SUNMemory buffer_data = NULL;
    SUNCheckCall(SUNMemoryHelper_Alloc(IMPL_PROP(self, mem_helper), &buffer_data,
                                       buffer_size, buffer_mem_type, queue));

    SUNCheckCall(SUNMemoryHelper_Copy(IMPL_PROP(self, mem_helper), buffer_data,
                                      leaf_data, buffer_size, queue));

    sunrealtype* data_ptr = leaf_data->ptr;
    *t                    = data_ptr[0];
    data_ptr              = buffer_data->ptr;

    SUNCheckCall(N_VBufUnpack(v, &data_ptr[1]));

    SUNCheckCall(
      SUNMemoryHelper_Dealloc(IMPL_PROP(self, mem_helper), buffer_data, queue));
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetData_InMem(SUNDataNode self, SUNMemoryType src_mem_type,
                                     SUNMemoryType node_mem_type, void* data,
                                     size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  void* queue = NULL;

  SUNAssert(BASE_PROP(self, dtype) == SUNDATANODE_LEAF, SUN_ERR_ARG_WRONGTYPE);

  SUNMemory data_mem_src = SUNMemoryHelper_Wrap(IMPL_PROP(self, mem_helper),
                                                data, src_mem_type);
  SUNCheckLastErr();

  SUNMemory data_mem_dst = NULL;
  SUNCheckCall(SUNMemoryHelper_AllocStrided(IMPL_PROP(self, mem_helper),
                                            &data_mem_dst, data_bytes,
                                            data_stride, node_mem_type, queue));

  SUNCheckCall(SUNMemoryHelper_Copy(IMPL_PROP(self, mem_helper), data_mem_dst,
                                    data_mem_src, data_bytes, queue));

  SUNMemoryHelper_Dealloc(IMPL_PROP(self, mem_helper), data_mem_src, queue);

  IMPL_PROP(self, leaf_data) = data_mem_dst;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetDataNvector_InMem(SUNDataNode self, N_Vector v,
                                            sunrealtype t)
{
  SUNFunctionBegin(self->sunctx);

  void* queue = NULL;

  SUNMemoryType leaf_mem_type = SUNMEMTYPE_HOST;
  SUNMemoryType buffer_mem_type = N_VGetDeviceArrayPointer(v) ? SUNMEMTYPE_DEVICE
                                                              : SUNMEMTYPE_HOST;

  sunindextype buffer_size = 0;
  SUNCheckCall(N_VBufSize(v, &buffer_size));

  /* We allocate 1 extra sunrealtype for storing t */
  SUNMemory leaf_data = NULL;
  SUNCheckCall(
    SUNMemoryHelper_AllocStrided(IMPL_PROP(self, mem_helper), &leaf_data,
                                 buffer_size + sizeof(sunrealtype),
                                 sizeof(sunrealtype), leaf_mem_type, queue));

  if (leaf_mem_type == buffer_mem_type)
  {
    sunrealtype* data_ptr = leaf_data->ptr;
    data_ptr[0]           = t;
    SUNCheckCall(N_VBufPack(v, &data_ptr[1]));
  }
  else
  {
    /* If the node memory type is not the same as the N_Vector's memory type,
       then we will first need to create a buffer of the same type as the N_Vector's
       and then copy it to the node data. */
    SUNMemory buffer_data = NULL;
    SUNCheckCall(
      SUNMemoryHelper_AllocStrided(IMPL_PROP(self, mem_helper), &buffer_data,
                                   buffer_size + sizeof(sunrealtype),
                                   sizeof(sunrealtype), buffer_mem_type, queue));

    sunrealtype* data_ptr = buffer_data->ptr;
    data_ptr[0]           = t;
    SUNCheckCall(N_VBufPack(v, &data_ptr[1]));

    SUNCheckCall(SUNMemoryHelper_Copy(IMPL_PROP(self, mem_helper), leaf_data,
                                      buffer_data,
                                      buffer_size + sizeof(sunrealtype), queue));

    SUNCheckCall(
      SUNMemoryHelper_Dealloc(IMPL_PROP(self, mem_helper), buffer_data, queue));
  }

  IMPL_PROP(self, leaf_data) = leaf_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_Destroy_InMem(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  void* queue = NULL;

  if (BASE_PROP(*node, dtype) == SUNDATANODE_OBJECT)
  {
    SUNHashMap map = IMPL_PROP(*node, named_children);
    SUNHashMap_Destroy(&map);
  }
  else if (BASE_PROP(*node, dtype) == SUNDATANODE_LIST)
  {
    SUNStlVector_SUNDataNode_Destroy(&IMPL_PROP(*node, anon_children));
  }
  else if (BASE_PROP(*node, dtype) == SUNDATANODE_LEAF)
  {
    if (IMPL_PROP(*node, leaf_data))
    {
      SUNCheckCall(SUNMemoryHelper_Dealloc(IMPL_PROP(*node, mem_helper),
                                           IMPL_PROP(*node, leaf_data), queue));
    }
  }

  sunDataNodeInMem_DestroyEmpty(node);
  *node = NULL;

  return SUN_SUCCESS;
}

/* This function is the callback provided to the child hashmap as the destroy function. */
static void sunDataNodeFreeKeyValue(SUNDIALS_MAYBE_UNUSED SUNHashMapKeyValue* kv_ptr)
{
  // if (!kv_ptr || !(*kv_ptr)) { return; }
  // SUNDataNode value = (SUNDataNode)(*kv_ptr)->value;
  // SUNDataNode_Destroy_InMem(&value);
  /* Do nothing. We want the user of the class to have to call SUNDataNode_Destroy
     for each SUNDataNode, even child nodes.*/
  return;
}

/* This function is the callback provided to the child stlvector as the destroy function. */
static void sunDataNodeFreeValue(SUNDIALS_MAYBE_UNUSED SUNDataNode* nodeptr)
{
  // if (!nodeptr || !(*nodeptr)) { return; }
  // SUNDataNode_Destroy_InMem(nodeptr);
  /* Do nothing. We want the user of the class to have to call SUNDataNode_Destroy
     for each SUNDataNode, even child nodes.*/
  return;
}
