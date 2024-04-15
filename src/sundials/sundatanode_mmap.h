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
#include "sundials_datanode.h"
#include "sundials_hashmap_impl.h"

#ifndef SUNDATANODE_MMAP_H_
#define SUNDATANODE_MMAP_H_

#ifdef __cplusplus
extern "C" {
#endif

#define TTYPE SUNDataNode
#include "sundials_arraylist.h"

typedef struct SUNDataNode_MmapImpl_s* SUNDataNode_MmapImpl;

struct SUNDataNode_MmapImpl_s {
  SUNDataNode parent;

  // Node can only be an object, leaf, or list. It cannot be more than one of these at a time.

  // Properties for Leaf nodes (nodes that store data)
  void* leaf_data;
  size_t data_stride;
  size_t data_bytes;

  // Properties for Object nodes (nodes that are a collection of named nodes)
  const char* name;
  SUNHashMap named_children;
  sundataindex_t num_named_children;
  sundataindex_t max_named_children;

  // Properties for a List node (nodes that are a collection of anonymous nodes)
  // SUNArrayList_SUNDataNode anon_children;
  SUNDataNode* anon_children;
  sundataindex_t num_anon_children;
  sundataindex_t max_anon_children;
};

SUNErrCode SUNDataNode_CreateList_Mmap(sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out);

SUNErrCode SUNDataNode_CreateObject_Mmap(sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node_out);

SUNErrCode SUNDataNode_CreateLeaf_Mmap(void* leaf_data, size_t data_stride, size_t data_bytes,
                                       SUNContext sunctx, SUNDataNode* node_out);

SUNErrCode SUNDataNode_IsLeaf_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_IsList_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_IsObject_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_HasChildren_Mmap(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNErrCode SUNDataNode_AddChild_Mmap(SUNDataNode parent_node, SUNDataNode child_node);

SUNErrCode SUNDataNode_AddNamedChild_Mmap(SUNDataNode parent_node, const char* name, SUNDataNode child_node);

SUNErrCode SUNDataNode_GetChild_Mmap(const SUNDataNode node, sundataindex_t index, SUNDataNode* child_node);

SUNErrCode SUNDataNode_GetNamedChild_Mmap(const SUNDataNode node, const char* name, SUNDataNode* child_node);

SUNErrCode SUNDataNode_RemoveChild_Mmap(SUNDataNode node, sundataindex_t index, SUNDataNode* child_node);

SUNErrCode SUNDataNode_RemoveNamedChild_Mmap(SUNDataNode node, const char* name, SUNDataNode* child_node);

SUNErrCode SUNDataNode_GetData_Mmap(const SUNDataNode node, void** data);

SUNErrCode SUNDataNode_SetData_Mmap(SUNDataNode node, void* data, size_t data_stride, size_t data_bytes);

SUNErrCode SUNDataNode_Destroy_Mmap(SUNDataNode* node);

#ifdef __cplusplus
}
#endif

#endif // SUNDATANODE_MMAP_H_
