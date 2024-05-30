/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * SUNDIALS macros
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_MACROS_H
#define _SUNDIALS_MACROS_H

#include "sundials/sundials_config.h"

/* -----------------------------------------------------------------------------
 * SUNDIALS_MAYBE_UNUSED
 *
 * This maps to an attribute that can be used to silence warnings about unused
 * classes, typedefs, variables, functions, or methods when the entity cannot be
 * removed. For example, functions or variables that are only used when error
 * checks or profiling is enabled.
 * ---------------------------------------------------------------------------*/

#if __cplusplus >= 201703L || __STDC_VERSION__ > 201710L
#define SUNDIALS_MAYBE_UNUSED [[maybe_unused]]
#elif defined(SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_UNUSED)
#define SUNDIALS_MAYBE_UNUSED __attribute__((unused))
#else
#define SUNDIALS_MAYBE_UNUSED
#endif

#endif /* _SUNDIALS_MACROS_H */
