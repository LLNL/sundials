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
 # To use the StlVector, first define TTYPE with your data type
 # before including this header. If you need StlVectors that hold
 # different types in the same file, then define TTYPE for the first,
 # include this header, then repeat.
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/sundials_core.h>

#ifndef TTYPE
#error "Must define template type for SUNStlVector"
#endif

#define CONCAT(a, b)            a##b
#define PASTE(a, b)             CONCAT(a, b)
#define MAKE_NAME(prefix, name) PASTE(prefix, PASTE(_, name))

#define SUNStlVectorTtype_s MAKE_NAME(SUNStlVector, PASTE(TTYPE, _s))
#define SUNStlVectorTtype   MAKE_NAME(SUNStlVector, TTYPE)

typedef struct SUNStlVectorTtype_s* SUNStlVectorTtype;

struct SUNStlVectorTtype_s
{
  int64_t size;
  int64_t capacity;
  TTYPE* values;
  void (*destroyValue)(TTYPE*);
};

// This constant controls how much space will be allocated when a resize is needed.
// The new capacity is GROWTH_FACTOR*current_capacity.
// Some std::vector implementations use 2, but 1.5 will be more conservative in terms
// of the memory usage but yields a larger constant factor in terms of the
// amortized constant time complexity.
#define GROWTH_FACTOR 1.5

static inline SUNStlVectorTtype MAKE_NAME(SUNStlVectorTtype,
                                          New)(int64_t init_capacity,
                                               void (*destroyValue)(TTYPE*))
{
  SUNStlVectorTtype self =
    (SUNStlVectorTtype)malloc(sizeof(struct SUNStlVectorTtype_s));
  self->size         = 0;
  self->capacity     = init_capacity > 0 ? init_capacity : 1;
  self->values       = (TTYPE*)malloc(sizeof(TTYPE) * self->capacity);
  self->destroyValue = destroyValue;
  return self;
}

static inline sunbooleantype MAKE_NAME(SUNStlVectorTtype,
                                       IsEmpty)(SUNStlVectorTtype self)
{
  return self->size == 0;
}

static inline void MAKE_NAME(SUNStlVectorTtype, Resize)(SUNStlVectorTtype self,
                                                        int64_t new_capacity)
{
  if (new_capacity <= self->capacity) return;
  TTYPE* new_values = (TTYPE*)realloc(self->values, sizeof(TTYPE) * new_capacity);
  self->values   = new_values;
  self->capacity = new_capacity;
}

static inline void MAKE_NAME(SUNStlVectorTtype, Grow)(SUNStlVectorTtype self)
{
  if (self->size == self->capacity)
  {
    int64_t new_capacity =
      (int64_t)(ceil(((double)self->capacity) * GROWTH_FACTOR));
    MAKE_NAME(SUNStlVectorTtype, Resize)(self, new_capacity);
  }
}

static inline void MAKE_NAME(SUNStlVectorTtype,
                             PushBack)(SUNStlVectorTtype self, TTYPE element)
{
  if (self->size == self->capacity)
  {
    MAKE_NAME(SUNStlVectorTtype, Grow)(self);
  }
  self->values[self->size++] = element;
}

static inline TTYPE* MAKE_NAME(SUNStlVectorTtype, At)(SUNStlVectorTtype self,
                                                      int64_t index)
{
  if (index >= self->size)
  {
    // Handle index out of bounds
    return NULL;
  }
  return &(self->values[index]);
}

static inline void MAKE_NAME(SUNStlVectorTtype, Set)(SUNStlVectorTtype self,
                                                     int64_t index, TTYPE element)
{
  if (index >= self->size)
  {
    // Handle index out of bounds
    return;
  }
  self->values[index] = element;
}

static inline void MAKE_NAME(SUNStlVectorTtype, PopBack)(SUNStlVectorTtype self)
{
  static TTYPE nullish;
  if (self->size == 0) return;
  self->size--;
  MAKE_NAME(SUNStlVectorTtype, Set)(self, self->size, nullish);
}

static inline void MAKE_NAME(SUNStlVectorTtype, Erase)(SUNStlVectorTtype self,
                                                       int64_t index)
{
  static TTYPE nullish;
  if (self->size == 0) return;
  MAKE_NAME(SUNStlVectorTtype, Set)(self, index, nullish);
  for (int64_t i = index; i < self->size - 1; i++)
  {
    self->values[i]     = self->values[i + 1];
    self->values[i + 1] = nullish;
  }
  self->size -= 1;
}

static inline int64_t MAKE_NAME(SUNStlVectorTtype, Size)(SUNStlVectorTtype self)
{
  return self->size;
}

static inline int64_t MAKE_NAME(SUNStlVectorTtype,
                                Capacity)(SUNStlVectorTtype self)
{
  return self->capacity;
}

static inline void MAKE_NAME(SUNStlVectorTtype,
                             Destroy)(SUNStlVectorTtype* self_ptr)
{
  static TTYPE nullish;

  if (!self_ptr || !(*self_ptr)) return;

  SUNStlVectorTtype self = *self_ptr;

  for (int64_t i = 0; i < MAKE_NAME(SUNStlVectorTtype, Size)(self); i++)
  {
    self->destroyValue(&(self->values[i]));
    self->values[i] = nullish;
  }

  *self_ptr = NULL;
}

#undef TTYPE
