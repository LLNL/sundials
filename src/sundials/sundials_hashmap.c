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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"

static const uint64_t HASH_PRIME        = 14695981039346656037U;
static const uint64_t HASH_OFFSET_BASIS = 1099511628211U;

/*
  For a nice discussion on popular hashing algorithms see:
  https://softwareengineering.stackexchange.com/questions/49550/which-hashing-algorithm-is-best-for-uniqueness-and-speed/145633#145633

  This is a 64-bit implementation of the 'a' modification of the
  Fowler-Noll-Vo hash (i.e., FNV1-a).
 */
static uint64_t fnv1a_hash(const char* str)
{
  uint64_t hash = HASH_OFFSET_BASIS;
  char c;
  while ((c = *str++)) { hash = (hash ^ c) * HASH_PRIME; }
  return hash;
}

/*
  This function creates a new SUNHashMap object allocated to hold
  up to 'capacity' entries.

  **Arguments:**
    * ``capacity`` -- the initial capactity number of the hashmap
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 allocated

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_New(size_t capacity, SUNHashMap* map)
{
  if (capacity <= 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  *map = NULL;
  *map = (SUNHashMap)malloc(sizeof(**map));

  if (!map) { return SUN_ERR_MALLOC_FAIL; }

  (*map)->capacity = capacity;

  SUNArrayList_SUNHashMapKeyValue buckets =
    SUNArrayList_SUNHashMapKeyValue_New(capacity);
  if (!buckets)
  {
    free(*map);
    return SUN_ERR_MALLOC_FAIL;
  }

  /* Initialize all buckets to NULL */
  for (size_t i = 0; i < capacity; i++)
  {
    SUNArrayList_SUNHashMapKeyValue_PushBack(buckets, NULL);
  }

  (*map)->buckets = buckets;

  return SUN_SUCCESS;
}

size_t SUNHashMap_Capacity(SUNHashMap map)
{
  return SUNArrayList_SUNHashMapKeyValue_Capacity(map->buckets);
}

/*
  This function frees the SUNHashMap object.

  **Arguments:**
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 deallocated and set to ``NULL``
    * ``freevalue`` -- callback function that should free the value object

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_Destroy(SUNHashMap* map, void (*freevalue)(void* ptr))
{
  if (map == NULL || freevalue == NULL) { return SUN_SUCCESS; }

  SUNArrayList_SUNHashMapKeyValue buckets = (*map)->buckets;
  for (size_t i = 0; i < (*map)->capacity; i++)
  {
    SUNHashMapKeyValue bucket = *SUNArrayList_SUNHashMapKeyValue_At(buckets, i);
    if (bucket && bucket->value) { freevalue(bucket); }
    if (bucket)
    {
      free(bucket->key);
      free(bucket);
    }
  }

  free(buckets);
  free(*map);
  *map = NULL;

  return SUN_SUCCESS;
}

/*
  This function iterates the map over the range [start, N]. N is either the
  index at which ``yieldfn`` indicates the iteration should stop, or the max
  entries in the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``start`` -- the start of the iteration range
    * ``yieldfn`` -- the callback function to call every iteration
                     this should return -1 to continue the iteration, or >= 0 to
                     stop; the first argument is the current index, the second
                     argument is the current key-value pair, and the final
                     argument is the same pointer ``ctx`` as the final argument
                     to SUNHashMapIterate.
    * ``ctx`` -- a pointer to pass on to ``yieldfn``

  **Returns:**
    * ``capacity + 1`` -- iterated the whole map
    * ``>=0`` -- the index at which the iteration stopped
    * ``<-1`` -- an error occurred
 */
size_t SUNHashMap_Iterate(SUNHashMap map, int start,
                          int (*yieldfn)(int, SUNHashMapKeyValue, void*),
                          void* ctx)
{
  if (map == NULL || yieldfn == NULL) { return (-2); }

  for (size_t i = start; i < SUNArrayList_SUNHashMapKeyValue_Size(map->buckets);
       i++)
  {
    int retval =
      yieldfn(i, *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, i), ctx);
    if (retval >= 0)
    {
      return (retval); /* yieldfn indicates the loop should break */
    }
    if (retval < -1) { return (retval); /* error occurred */ }
  }

  return SUNArrayList_SUNHashMapKeyValue_Size(map->buckets) + 1;
}

static int sunHashMapLinearProbeInsert(int idx, SUNHashMapKeyValue kv,
                                       SUNDIALS_MAYBE_UNUSED const void* ctx)
{
  /* find the next open spot */
  if (kv == NULL) { return (idx); /* open spot found at idx */ }
  return (-1); /* keep looking */
}

/*
  This function creates a key-value pair and attempts to insert it into the map.
  Will use linear probing if there is a collision.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to store (we will make a copy)
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
    * ``-2`` -- duplicate key
 */
int SUNHashMap_Insert(SUNHashMap map, const char* key, void* value)
{
  size_t idx;
  size_t retval;
  SUNHashMapKeyValue kvp;

  if (map == NULL || key == NULL || value == NULL) { return (-1); }

  /* We want the index to be in (0, SUNHashMap_Capacity(map)) */
  idx = (size_t)(fnv1a_hash(key) % SUNHashMap_Capacity(map));

  /* Check if the bucket is already filled (i.e., we might have had a collision) */
  kvp = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, idx);
  if (kvp != NULL)
  {
    /* Determine if key is actually a duplicate (not allowed) */
    if (!strcmp(key, kvp->key)) { return (-2); }

    /* OK, it was a real collision, so find the next open spot */
    retval = SUNHashMap_Iterate(map, idx, sunHashMapLinearProbeInsert, NULL);
    if (retval < 0) { return (-1); /* error occurred */ }
    if (retval >= SUNHashMap_Capacity(map))
    {
      /* There are no open spaces, so resize the hashmap and then try to insert again. */
      size_t old_capacity = SUNHashMap_Capacity(map);
      SUNArrayList_SUNHashMapKeyValue_Grow(map->buckets);
      /* Set all of the new possible elements in the ArrayList to NULL
         because we always want to have the list capacity == list size. */
      for (size_t i = old_capacity; i < SUNHashMap_Capacity(map); i++)
      {
        SUNArrayList_SUNHashMapKeyValue_PushBack(map->buckets, NULL);
      }
      return SUNHashMap_Insert(map, key, value);
    }

    idx = retval;
  }

  /* Create the key-value pair */
  kvp = (SUNHashMapKeyValue)malloc(sizeof(*kvp));

  /* Copy the original_key so that the hashmap owns it */
  size_t len     = strlen(key);
  char* key_copy = malloc(sizeof(*key) * len);
  strcpy(key_copy, key);

  kvp->key   = key_copy;
  kvp->value = value;

  /* Insert the key-value pair */
  SUNArrayList_SUNHashMapKeyValue_Set(map->buckets, idx, kvp);

  return (0);
}

static int sunHashMapLinearProbeGet(int idx, SUNHashMapKeyValue kv,
                                    const void* key)
{
  /* target key cannot be NULL */
  if (key == NULL) { return (-2); }

  /* find the matching entry */
  if (kv == NULL) { return (-1); /* keep looking since this bucket is empty */ }
  if (!strcmp(kv->key, (const char*)key))
  {
    return (idx); /* found it at idx */
  }
  return (-1); /* keep looking */
}

/*
  This function gets the value for the given key.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to look up
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
    * ``-2`` -- key not found
 */
int SUNHashMap_GetValue(SUNHashMap map, const char* key, void** value)
{
  size_t idx;
  size_t retval;

  if (map == NULL || key == NULL || value == NULL) { return (-1); }

  /* We want the index to be in (0, SUNHashMap_Capacity(map)) */
  idx = (int)(fnv1a_hash(key) % SUNHashMap_Capacity(map));

  SUNHashMapKeyValue kvp = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, idx);

  /* Check if the key exists */
  if (kvp == NULL) { return (-2); }

  /* Check to see if this is a collision */
  if (strcmp(kvp->key, key))
  {
    /* Keys did not match, so we have a collision and need to probe */
    retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeGet, key);
    if (retval < 0) { return (-1); /* error occurred */ }
    if (retval > SUNHashMap_Capacity(map)) { return (-2); /* not found */ }
    else { idx = retval; }
  }

  /* Return a reference to the value only */
  kvp    = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, idx);
  *value = kvp->value;

  return (0);
}

/*
  This function remove the key-value pair.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to remove
    * ``value`` -- the value to remove

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
    * ``-2`` -- key not found
 */
int SUNHashMap_Remove(SUNHashMap map, const char* key, void** value)
{
  size_t idx;
  size_t retval;

  if (map == NULL || key == NULL) { return (-1); }

  /* We want the index to be in (0, SUNHashMap_Capacity(map)) */
  idx = (int)(fnv1a_hash(key) % SUNHashMap_Capacity(map));

  SUNHashMapKeyValue kvp = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, idx);

  /* Check if the key exists */
  if (kvp == NULL) { return (-2); }

  /* Check to see if this is a collision */
  if (strcmp(kvp->key, key))
  {
    /* Keys did not match, so we have a collision and need to probe */
    retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeGet,
                                (void*)key);
    if (retval < 0) { return (-1); /* error occurred */ }
    else if (retval > SUNHashMap_Capacity(map)) { return (-2); /* not found */ }
    else { idx = retval; }
  }

  /* Return a reference to the value only */
  kvp    = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, idx);
  *value = kvp->value;

  /* Clear the bucket by setting it to NULL */
  SUNArrayList_SUNHashMapKeyValue_Set(map->buckets, idx, NULL);

  return (0);
}

/*
  This function allocates a new array the same max_size as the map,
  then it sorts map into a new array of key-value pairs leaving
  the map unchanged.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``sorted`` -- pointer to the sorted array of key-value pairs, this
                    function will allocate the array
    * ``compar`` -- comparator function that is passed to the C standard qsort
                    function

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_Sort(SUNHashMap map, SUNHashMapKeyValue** sorted,
                           int (*compar)(const void*, const void*))
{
  int i;

  if (!map || !compar) { return SUN_ERR_ARG_CORRUPT; }

  *sorted =
    (SUNHashMapKeyValue*)malloc(SUNHashMap_Capacity(map) * sizeof(**sorted));
  if (!(*sorted)) { return SUN_ERR_MALLOC_FAIL; }

  /* Copy the buckets into a new array */
  for (i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    (*sorted)[i] = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, i);
  }

  qsort(*sorted, SUNHashMap_Capacity(map), sizeof(SUNHashMapKeyValue), compar);

  return SUN_SUCCESS;
}

/*
  This function allocates a new array with just they values of the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``values`` -- pointer to the array of keys
    * ``value_size`` -- the size of the values in bytes

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_Values(SUNHashMap map, void*** values, size_t value_size)
{
  int i;
  int count = 0;

  if (!map) { return SUN_ERR_ARG_CORRUPT; }

  *values = (void**)malloc(SUNHashMap_Capacity(map) * value_size);
  if (!values) { return SUN_ERR_MALLOC_FAIL; }

  /* Copy the values into a new array */
  for (i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    SUNHashMapKeyValue kvp = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, i);
    if (kvp) { (*values)[count++] = kvp->value; }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNHashMap_PrintKeys(SUNHashMap map, FILE* file)
{
  int i;
  int count = 0;

  if (!map) { return SUN_ERR_ARG_CORRUPT; }

  /* Print keys into a new array */
  for (i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    SUNHashMapKeyValue kvp = *SUNArrayList_SUNHashMapKeyValue_At(map->buckets, i);
    if (kvp) { fprintf(file, "%s\n", kvp->key); }
  }

  return SUN_SUCCESS;
}
