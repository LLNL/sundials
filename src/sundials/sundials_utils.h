/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file contains common utility functions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_UTILS_H
#define _SUNDIALS_UTILS_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"

/* ----------------------------------------- *
 * Vector creation and destruction utilities *
 * ----------------------------------------- */

static inline SUNErrCode sunVec_Clone(N_Vector tmpl, N_Vector* v)
{
  if (*v != NULL) { return SUN_SUCCESS; }
  *v = N_VClone(tmpl);
  if (*v == NULL) { return SUN_ERR_MEM_FAIL; }
  tmpl->sunctx->vec_count++;
  return SUN_SUCCESS;
}

static inline SUNErrCode sunVec_Destroy(N_Vector* v)
{
  if (v == NULL) { return SUN_SUCCESS; }
  if (*v == NULL) { return SUN_SUCCESS; }
  SUNContext sunctx = (*v)->sunctx;
  N_VDestroy(*v);
  *v = NULL;
  sunctx->vec_count--;
  return SUN_SUCCESS;
}

static inline SUNErrCode sunVecArray_Clone(int count, N_Vector tmpl,
                                           N_Vector** v)
{
  if (count < 1) { return SUN_ERR_ARG_OUTOFRANGE; }
  if (*v != NULL) { return SUN_SUCCESS; }
  *v = N_VCloneVectorArray(count, tmpl);
  if (*v == NULL) { return SUN_ERR_MEM_FAIL; }
  tmpl->sunctx->vec_count += count;
  return SUN_SUCCESS;
}

static inline SUNErrCode sunVecArray_Destroy(int count, N_Vector** v)
{
  if (count < 1) { return SUN_ERR_ARG_OUTOFRANGE; }
  if (v == NULL) { return SUN_SUCCESS; }
  if (*v == NULL) { return SUN_SUCCESS; }
  SUNContext sunctx = ((*v)[0])->sunctx;
  N_VDestroyVectorArray(*v, count);
  *v = NULL;
  sunctx->vec_count -= count;
  return SUN_SUCCESS;
}

/* ------------------ *
 * Printing utilities *
 * ------------------ */

/* width of name field in sunfprintf_<type> for aligning table output */
#define SUN_TABLE_WIDTH 29

static inline char* sunSignedToString(int64_t val)
{
  char* str     = NULL;
  size_t length = snprintf(NULL, 0, "%lld", (long long)val);
  str           = (char*)malloc(sizeof(*str) * (length + 1));
  snprintf(str, length + 1, "%lld", (long long)val);
  return str;
}

static inline char* sunCombineFileAndLine(int line, const char* file)
{
  size_t total_str_len = strlen(file) + 6;
  char* file_and_line  = (char*)malloc(total_str_len * sizeof(char));
  snprintf(file_and_line, total_str_len, "%s:%d", file, line);
  return file_and_line;
}

/*
 * Implementation of the GNU extension function vasprintf which
 * is itself an analog for vsprintf, except it allocates a string
 * large enough to hold the output byte ('\0').
 */
static inline int sunvasnprintf(char** str, const char* fmt, va_list args)
{
  int size = 0;

  /* compute string length */
  va_list tmp1;
  va_copy(tmp1, args);
  size = vsnprintf(NULL, 0, fmt, tmp1);
  va_end(tmp1);

  if (size < 0) { return -1; }

  /* add one to size for the null terminator*/
  *str = (char*)malloc(size + 1);
  if (NULL == *str) { return -1; }

  va_list tmp2;
  va_copy(tmp2, args);
  size = vsnprintf(*str, size + 1, fmt, tmp2);
  va_end(tmp2);

  return size;
}

static inline void sunCompensatedSum(sunrealtype base, sunrealtype inc,
                                     sunrealtype* sum, sunrealtype* error)
{
  sunrealtype err           = *error;
  volatile sunrealtype tmp1 = inc - err;
  volatile sunrealtype tmp2 = base + tmp1;
  *error                    = (tmp2 - base) - tmp1;
  *sum                      = tmp2;
}

static inline void sunfprintf_real(FILE* fp, SUNOutputFormat fmt,
                                   sunbooleantype start, const char* name,
                                   sunrealtype value)
{
  if (fmt == SUN_OUTPUTFORMAT_TABLE)
  {
    fprintf(fp, "%-*s = " SUN_FORMAT_G "\n", SUN_TABLE_WIDTH, name, value);
  }
  else
  {
    if (!start) { fprintf(fp, ","); }
    fprintf(fp, "%s," SUN_FORMAT_E, name, value);
  }
}

static inline void sunfprintf_long(FILE* fp, SUNOutputFormat fmt,
                                   sunbooleantype start, const char* name,
                                   long value)
{
  if (fmt == SUN_OUTPUTFORMAT_TABLE)
  {
    fprintf(fp, "%-*s = %ld\n", SUN_TABLE_WIDTH, name, value);
  }
  else
  {
    if (!start) { fprintf(fp, ","); }
    fprintf(fp, "%s,%ld", name, value);
  }
}

static inline void sunfprintf_long_array(FILE* fp, SUNOutputFormat fmt,
                                         sunbooleantype start, const char* name,
                                         long* value, size_t count)
{
  if (count < 1) { return; }

  if (fmt == SUN_OUTPUTFORMAT_TABLE)
  {
    fprintf(fp, "%-*s = %ld", SUN_TABLE_WIDTH, name, value[0]);
    for (size_t i = 1; i < count; i++) { fprintf(fp, ", %ld", value[i]); }
    fprintf(fp, "\n");
  }
  else
  {
    if (!start) { fprintf(fp, ","); }
    for (size_t i = 0; i < count; i++)
    {
      fprintf(fp, "%s %zu,%ld", name, i, value[i]);
    }
  }
}

/* ------------ *
 * Vector stack *
 * ------------ */

struct SUNVecStack_
{
  int max;        /* maximum size of the stack */
  int top;        /* index of top of the stack, -1 = empty */
  N_Vector tmpl;  /* template vector to clone from */
  N_Vector* vecs; /* stack of temporary vectors */
};

static inline SUNErrCode SUNVecStack_Create(N_Vector tmpl,
                                            SUNVecStack* stack_out)
{
  /* SUNAssert(tmpl, SUN_ERR_ARG_CORRUPT); */
  /* SUNAssert(stack_out, SUN_ERR_ARG_CORRUPT); */

  SUNVecStack stack = (SUNVecStack)malloc(sizeof(struct SUNVecStack_));
  /* SUNAssert(stack, SUN_ERR_MALLOC_FAIL); */

  stack->max  = 0;
  stack->top  = -1;
  stack->tmpl = tmpl;
  stack->vecs = NULL;

  *stack_out = stack;

  return SUN_SUCCESS;
}

static inline SUNErrCode SUNVecStack_Destroy(SUNVecStack* stack_in)
{
  if (stack_in == NULL) { return SUN_SUCCESS; }
  if (*stack_in == NULL) { return SUN_SUCCESS; }

  /* SUNCheck(sunVecArray_Destroy((*stack_in)->max, &((*stack)->vecs)), */
  /*          SUN_ERR_DESTROY_FAIL); */
  sunVecArray_Destroy((*stack_in)->max, &((*stack_in)->vecs));
  free(*stack_in);
  *stack_in = NULL;

  return SUN_SUCCESS;
}

static inline SUNErrCode SUNVecStack_Pop(SUNVecStack stack, N_Vector* vec_out)
{
  /* SUNAssert(stack, SUN_ERR_ARG_CORRUPT); */
  /* SUNAssert(vec_out, SUN_ERR_ARG_CORRUPT); */
  /* SUNAssert(stack->top >= -1, SUN_ERR_OUTOFRANGE); */
  /* SUNAssert(stack->top < stack->max, SUN_ERR_OUTOFRANGE); */

  if (stack->top < 0)
  {
    /* create new stack of vectors */
    if (stack->vecs) { free(stack->vecs); }
    stack->vecs = (N_Vector*)malloc((stack->max + 1) * sizeof(N_Vector));
    /* SUNAssert(stack->vecs, SUN_ERR_MALLOC_FAIL); */
    stack->max++;

    /* create new vector on the bottom of the stack */
    /* SUNCheck(sunVec_Clone(sunctx, stack->tmpl, &(stack->vecs[0])), */
    /*          SUN_ERR_MEMFAIL); */
    sunVec_Clone(stack->tmpl, &(stack->vecs[0]));
    stack->top++;
  }

  /* remove vector from top of stack */
  *vec_out = stack->vecs[stack->top];
  stack->top--;

  return SUN_SUCCESS;
}

static inline SUNErrCode SUNVecStack_Push(SUNVecStack stack, N_Vector* vec_in)
{
  /* SUNAssert(stack, SUN_ERR_ARG_CORRUPT); */
  /* SUNAssert(vec_in, SUN_ERR_ARG_CORRUPT); */
  /* SUNAssert(stack->top >= -1, SUN_ERR_OUTOFRANGE); */
  /* SUNAssert(stack->top + 1 < stack->max, SUN_ERR_OUTOFRANGE); /\* full stack *\/ */

  /* add vector to the top of the stack */
  stack->top++;
  stack->vecs[stack->top] = *vec_in;

  /* user has released the vector, nullify input vector */
  *vec_in = NULL;

  return SUN_SUCCESS;
}

#endif /* _SUNDIALS_UTILS_H */
