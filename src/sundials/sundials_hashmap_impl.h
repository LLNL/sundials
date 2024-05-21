/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * A simple hashmap implementation for char* keys and
 * void* values. Uses linear probing to resolve collisions.
 * The values can be anything, but will be freed by
 * the hash map upon its destruction.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_HASHMAP_IMPL_H
#define _SUNDIALS_HASHMAP_IMPL_H

#include <stdlib.h>
#include <sundials/sundials_types.h>

typedef struct SUNHashMapKeyValue_* SUNHashMapKeyValue;

struct SUNHashMapKeyValue_
{
  const char* key;
  void* value;
};

typedef struct SUNHashMap_* SUNHashMap;

struct SUNHashMap_
{
  int size;     /* current number of entries */
  int max_size; /* max number of entries */
  SUNHashMapKeyValue* buckets;
};

SUNErrCode SUNHashMap_New(int max_size, SUNHashMap* map);
SUNErrCode SUNHashMap_Destroy(SUNHashMap* map, void (*freevalue)(void* ptr));
int SUNHashMap_Iterate(SUNHashMap map, int start,
                       int (*yieldfn)(int, SUNHashMapKeyValue, const void*),
                       const void* ctx);
int SUNHashMap_Insert(SUNHashMap map, const char* key, void* value);
int SUNHashMap_GetValue(SUNHashMap map, const char* key, void** value);
SUNErrCode SUNHashMap_Sort(SUNHashMap map, SUNHashMapKeyValue** sorted,
                           int (*compar)(const void*, const void*));

#if SUNDIALS_MPI_ENABLED
SUNErrCode SUNHashMap_Values(SUNHashMap map, void*** values, size_t value_size);
#endif

#endif
