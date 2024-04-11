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
 * Implementation of a resizable container similar to a std::vector.
 * The values can be anything but data must be contiguous.
 #
 # To use the ArrayList, first define TTYPE with your data type
 # before including this header. If you need ArrayLists that hold
 # different types in the same file, then define TTYPE for the first,
 # include this header, then repeat.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_core.h>
#include <stdlib.h>

#ifndef TTYPE
#error "Must define template type for SUNArrayList"
#endif

#define CONCAT(a, b) a##b
#define PASTE(a, b) CONCAT(a, b)
#define MAKE_NAME(prefix, name) PASTE(prefix, PASTE(_, name))

#define SUNArrayListTtype_s MAKE_NAME(SUNArrayList, PASTE(TTYPE, _s))
#define SUNArrayListTtype MAKE_NAME(SUNArrayList, TTYPE)

typedef struct SUNArrayListTtype_s* SUNArrayListTtype;

struct SUNArrayListTtype_s {
  size_t size;
  size_t capacity;
  TTYPE* values;
};

// This constant controls how much space will be allocated when a resize is needed.
// The new capacity is GROWTH_FACTOR*current_capacity.
// Some std::vector implementations use 2, but 1.5 will be more conservative in terms
// of the memory usage but yields a larger constant factor in terms of the
// amorirtized constant time complexity.
#define GROWTH_FACTOR 1.5

static inline
SUNArrayListTtype MAKE_NAME(SUNArrayListTtype, New)(size_t init_capacity) {
  SUNArrayListTtype list = (SUNArrayListTtype) malloc(sizeof(struct SUNArrayListTtype_s));
  list->size = 0;
  list->capacity = init_capacity > 0 ? init_capacity : 1;
  list->values = (TTYPE*) malloc(sizeof(TTYPE) * list->capacity);
  return list;
}

static inline
void MAKE_NAME(SUNArrayListTtype, Destroy)(SUNArrayListTtype* list) {
  if (!(*list)) return;
  free((*list)->values);
  free(*list);
  *list = NULL;
}

static inline
sunbooleantype MAKE_NAME(SUNArrayListTtype, IsEmpty)(SUNArrayListTtype list) {
  return list->size == 0;
}

static inline
void MAKE_NAME(SUNArrayListTtype, Resize)(SUNArrayListTtype list, size_t new_capacity) {
  if (new_capacity <= list->capacity) return;
  TTYPE* new_values = (TTYPE*) realloc(list->values, sizeof(TTYPE) * new_capacity);
  list->values = new_values;
  list->capacity = new_capacity;
}

static inline
void MAKE_NAME(SUNArrayListTtype, Grow)(SUNArrayListTtype list) {
  if (list->size == list->capacity) {
    size_t new_capacity = (size_t)(ceil(list->capacity * GROWTH_FACTOR));
    MAKE_NAME(SUNArrayListTtype, Resize)(list, new_capacity);
  }
}

static inline
void MAKE_NAME(SUNArrayListTtype, PushBack)(SUNArrayListTtype list, TTYPE element) {
  if (list->size == list->capacity) {
    MAKE_NAME(SUNArrayListTtype, Grow)(list);
  }
  list->values[list->size++] = element;
}

static inline
TTYPE* MAKE_NAME(SUNArrayListTtype, At)(SUNArrayListTtype list, int index) {
  if (index < 0 || index >= list->size) {
    // Handle index out of bounds
    return NULL;
  }
  return &(list->values[index]);
}

static inline
void MAKE_NAME(SUNArrayListTtype, Set)(SUNArrayListTtype list, int index, TTYPE element) {
  if (index < 0 || index >= list->size) {
    // Handle index out of bounds
    return;
  }
  list->values[index] = element;
}

static inline
void MAKE_NAME(SUNArrayListTtype, PopBack)(SUNArrayListTtype list) {
  if (list->size == 0) return;
  list->size--;
}

static inline
size_t MAKE_NAME(SUNArrayListTtype, Size)(SUNArrayListTtype list) {
  return list->size;
}

static inline
size_t MAKE_NAME(SUNArrayListTtype, Capacity)(SUNArrayListTtype list) {
  return list->capacity;
}

#undef TTYPE
