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
  SUNDATANODE_LEAF,
  SUNDATANODE_LIST,
  SUNDATANODE_OBJECT
} SUNDataNodeType;

typedef enum {
  SUNDATAIOMODE_MMAP,
} SUNDataIOMode;

typedef struct SUNDataNode_s* SUNDataNode;

struct SUNDataNode_s {
  SUNErrCode (*hasChildren)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*isLeaf)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*isList)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*isObject)(const SUNDataNode, sunbooleantype* yes_or_no);
  SUNErrCode (*addChild)(SUNDataNode, SUNDataNode child_node);
  SUNErrCode (*addNamedChild)(SUNDataNode, const char* name, SUNDataNode child_node);
  SUNErrCode (*getChild)(const SUNDataNode, sundataindex_t index, SUNDataNode* child_node);
  SUNErrCode (*getNamedChild)(const SUNDataNode, const char* name, SUNDataNode* child_node);
  SUNErrCode (*removeChild)(SUNDataNode, sundataindex_t index, SUNDataNode* child_node);
  SUNErrCode (*getData)(const SUNDataNode, void** data);
  SUNErrCode (*setData)(SUNDataNode, void* data, size_t data_stride, size_t data_bytes);
  SUNErrCode (*destroy)(SUNDataNode*);

  void* impl;
  SUNDataNodeType dtype;
  SUNContext sunctx;
};

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateEmpty(SUNContext sunctx, SUNDataNode* node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateLeaf(SUNDataIOMode io_mode, void* leaf_data, size_t data_stride, size_t data_bytes,
                                  SUNContext sunctx, SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateList(SUNDataIOMode io_mode, sundataindex_t num_elements,
                                  SUNContext sunctx, SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_CreateObject(SUNDataIOMode io_mode, sundataindex_t num_elements,
                                    SUNContext sunctx, SUNDataNode* node_out);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_IsLeaf(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_IsList(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_HasChildren(const SUNDataNode node, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_AddChild(SUNDataNode parent_node, SUNDataNode child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_AddNamedChild(SUNDataNode parent_node, const char* name, SUNDataNode child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetChild(const SUNDataNode parent_node, sundataindex_t index, SUNDataNode* child_node);

SUNDIALS_EXPORT
SUNErrCode SUNDataNode_GetNamedChild(const SUNDataNode parent_node, const char* name, SUNDataNode* child_node);

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
