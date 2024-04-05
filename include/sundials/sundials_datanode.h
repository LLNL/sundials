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
 * -----------------------------------------------------------------
 * SUNDIALS data store class. The data store manages a collection
 * of checkpoints saved to some storage device.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_DATANODE_H
#define _SUNDIALS_DATANODE_H

#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef sunindextype sundataindex_t;

typedef enum {
  SUNDATANODE_BASIC,
} SUNDataNode_ID;

typedef enum {
  SUNIOPROTOCOL_MMAP,
} SUNIOProtocol;

typedef struct SUNDataNode_s* SUNDataNode;

struct SUNDataNode_s {
  SUNErrCode (*hasChildren)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*addChild)(SUNDataNode, sundataindex_t index);
  SUNErrCode (*getChild)(const SUNDataNode, sundataindex_t index, SUNDataNode*);
  SUNErrCode (*removeChild)(SUNDataNode, sundataindex_t index);
  SUNErrCode (*getData)(const SUNDataNode, void** data);
  SUNErrCode (*setData)(SUNDataNode, void* data, size_t data_stride, size_t data_bytes);
  SUNErrCode (*destroy)(SUNDataNode*);

  // Node basics
  SUNDataNode parent;

  // Node can only be an object, leaf, or list. It cannot be more than one of these at a time.

  // Properties for Leaf nodes (nodes that store data)
  void* leaf_data;
  size_t data_stride;
  size_t data_bytes;

  // // Properties for Object nodes (nodes that are a collection of named nodes)
  // const char* name;
  // SUNHashMap* named_children;
  // sundataindex_t num_named_children;

  // Properties for a List node (nodes that are a collection of anonymous nodes)
  SUNDataNode* anon_children;
  sundataindex_t num_anon_children;
  sundataindex_t max_anon_children;

  void* impl;
  SUNContext sunctx;
};

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateList(sundataindex_t num_elements, SUNContext sunctx, SUNDataNode* node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateLeaf(void* leaf_data, size_t data_stride, size_t data_bytes, SUNContext sunctx, SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_NodeIsLeaf(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_NodeIsList(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_HasChildren(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_AddChild(SUNDataNode parent_node, SUNDataNode child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetChild(const SUNDataNode parent_node, sundataindex_t index, SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_RemoveChild(SUNDataNode node, sundataindex_t index, SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetData(const SUNDataNode node, void** data);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_SetData(SUNDataNode node, void* data, size_t data_stride, size_t data_bytes);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_Destroy(SUNDataNode* node);


#ifdef __cplusplus
}
#endif

#endif
