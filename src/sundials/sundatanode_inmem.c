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

#include <sundials/sundatanode_inmem.h>

#include "sundials/sundials_errors.h"
#include "sundials_hashmap_impl.h"

#define GET_IMPL(node)        ((SUNDataNode_InMemImpl)(node)->impl)
#define IMPL_PROP(node, prop) (GET_IMPL(node)->prop)
#define BASE_PROP(node, prop) ((node)->prop)

static SUNDataNode sunDataNodeMmap_CreateEmpty(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node;
  SUNCheckCallNoRet(SUNDataNode_CreateEmpty(sunctx, &node));

  node->hasChildren      = SUNDataNode_HasChildren_InMem;
  node->isLeaf           = SUNDataNode_IsLeaf_InMem;
  node->isList           = SUNDataNode_IsList_InMem;
  node->isObject         = SUNDataNode_IsObject_InMem;
  node->addChild         = SUNDataNode_AddChild_InMem;
  node->addNamedChild    = SUNDataNode_AddNamedChild_InMem;
  node->getChild         = SUNDataNode_GetChild_InMem;
  node->getNamedChild    = SUNDataNode_GetNamedChild_InMem;
  node->removeChild      = SUNDataNode_RemoveChild_InMem;
  node->removeNamedChild = SUNDataNode_RemoveNamedChild_InMem;
  node->getData          = SUNDataNode_GetData_InMem;
  node->setData          = SUNDataNode_SetData_InMem;
  node->destroy          = SUNDataNode_Destroy_InMem;

  SUNDataNode_InMemImpl impl =
    (SUNDataNode_InMemImpl)malloc(sizeof(struct SUNDataNode_InMemImpl_s));
  SUNAssertNoRet(impl, SUN_ERR_MEM_FAIL);

  impl->parent             = NULL;
  impl->leaf_data          = NULL;
  impl->data_stride        = 0;
  impl->data_bytes         = 0;
  impl->name               = NULL;
  impl->named_children     = NULL;
  impl->num_named_children = 0;
  impl->anon_children      = NULL;

  node->impl = (void*)impl;

  return node;
}

static void sunDataNodeMmap_DestroyEmpty(SUNDataNode* node)
{
  if (!node || !(*node)) { return; }
  if (BASE_PROP(*node, impl)) { free(BASE_PROP(*node, impl)); }
  BASE_PROP(*node, impl) = NULL;
}

SUNErrCode SUNDataNode_CreateList_InMem(sundataindex_t init_size,
                                        SUNContext sunctx, SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeMmap_CreateEmpty(sunctx);

  BASE_PROP(node, dtype)         = SUNDATANODE_LIST;
  IMPL_PROP(node, anon_children) = SUNArrayList_SUNDataNode_New(init_size);

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateObject_InMem(sundataindex_t init_size,
                                          SUNContext sunctx,
                                          SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeMmap_CreateEmpty(sunctx);

  BASE_PROP(node, dtype) = SUNDATANODE_OBJECT;

  SUNHashMap map;
  SUNCheckCall(SUNHashMap_New(init_size, &map));

  IMPL_PROP(node, named_children) = map;

  *node_out = node;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_CreateLeaf_InMem(void* leaf_data, size_t data_stride,
                                        size_t data_bytes, SUNContext sunctx,
                                        SUNDataNode* node_out)
{
  SUNFunctionBegin(sunctx);

  SUNDataNode node = sunDataNodeMmap_CreateEmpty(sunctx);

  BASE_PROP(node, dtype) = SUNDATANODE_LEAF;

  IMPL_PROP(node, leaf_data)   = leaf_data;
  IMPL_PROP(node, data_stride) = data_stride;
  IMPL_PROP(node, data_bytes)  = data_bytes;

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
  if (IMPL_PROP(self, anon_children)) {}
  *yes_or_no =
    (IMPL_PROP(self, anon_children) &&
     SUNArrayList_SUNDataNode_Size(IMPL_PROP(self, anon_children)) != 0) ||
    IMPL_PROP(self, num_named_children) != 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_AddChild_InMem(SUNDataNode self, SUNDataNode child_node)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_PROP(self, dtype) == SUNDATANODE_LIST, SUN_ERR_ARG_WRONGTYPE);

  SUNArrayList_SUNDataNode_PushBack(IMPL_PROP(self, anon_children), child_node);
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
    return SUN_ERR_DATANODE_MAXCHILDREN;
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

  SUNArrayList_SUNDataNode children = IMPL_PROP(self, anon_children);
  sunbooleantype has_children;
  SUNCheckCall(SUNDataNode_HasChildren_InMem(self, &has_children));

  *child_node = NULL;
  if (has_children)
  {
    *child_node = *SUNArrayList_SUNDataNode_At(children, index);
  }

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
      return SUN_ERR_DATANODE_NODENOTFOUND;
    }
  }
  else { return SUN_ERR_ARG_INCOMPATIBLE; }

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
    SUNArrayList_SUNDataNode_At(IMPL_PROP(self, anon_children), index);
  if (child_node_ptr)
  {
    *child_node = *child_node_ptr;
    if (*child_node)
    {
      IMPL_PROP(*child_node, parent) = NULL;
      SUNArrayList_SUNDataNode_Erase(IMPL_PROP(self, anon_children), index);
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

SUNErrCode SUNDataNode_GetData_InMem(const SUNDataNode self, void** data)
{
  SUNFunctionBegin(self->sunctx);

  *data = IMPL_PROP(self, leaf_data);

  return SUN_SUCCESS;
}

SUNErrCode SUNDataNode_SetData_InMem(SUNDataNode self, void* data,
                                     size_t data_stride, size_t data_bytes)
{
  SUNFunctionBegin(self->sunctx);

  SUNAssert(BASE_PROP(self, dtype) == SUNDATANODE_LEAF, SUN_ERR_ARG_WRONGTYPE);

  IMPL_PROP(self, leaf_data)   = data;
  IMPL_PROP(self, data_stride) = data_stride;
  IMPL_PROP(self, data_bytes)  = data_bytes;

  return SUN_SUCCESS;
}

static void sunHashMapFreeDataNode(void* nodeptr)
{
  SUNDataNode node = (SUNDataNode)nodeptr;
  SUNDataNode_Destroy_InMem(&node);
}

SUNErrCode SUNDataNode_Destroy_InMem(SUNDataNode* node)
{
  SUNFunctionBegin((*node)->sunctx);

  if (BASE_PROP(*node, dtype) == SUNDATANODE_OBJECT)
  {
    SUNHashMap map = IMPL_PROP(*node, named_children);
    SUNHashMap_Destroy(&map, sunHashMapFreeDataNode);
  }
  else if (BASE_PROP(*node, dtype) == SUNDATANODE_LIST)
  {
    free(IMPL_PROP(*node, anon_children));
  }
  // sunDataNodeMmap_DestroyEmpty(node);

  return SUN_SUCCESS;
}
