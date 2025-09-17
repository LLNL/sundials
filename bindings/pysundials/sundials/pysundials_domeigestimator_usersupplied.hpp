/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#ifndef _PYSUNDIALS_SUNDOMEIGESTIMATOR_USERSUPPLIED_HPP
#define _PYSUNDIALS_SUNDOMEIGESTIMATOR_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

#include <sundials/sundials_domeigestimator.hpp>

#include "pysundials_helpers.hpp"

namespace nb = nanobind;

using namespace sundials::experimental;

struct SUNDomEigEstimatorFunctionTable
{
  nb::object atimes;
};

inline SUNDomEigEstimatorFunctionTable* SUNDomEigEstimatorFunctionTable_Alloc()
{
  // We must use malloc since ARKodeFree calls free
  auto fn_table = static_cast<SUNDomEigEstimatorFunctionTable*>(
    std::malloc(sizeof(SUNDomEigEstimatorFunctionTable)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(SUNDomEigEstimatorFunctionTable));

  return fn_table;
}

SUNErrCode sundomeigestimator_atimes_wrapper(void* A_data, N_Vector v, N_Vector z)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<SUNATimesFn>, SUNDomEigEstimatorFunctionTable,
    3>(&SUNDomEigEstimatorFunctionTable::atimes, A_data, v, z);
}

#endif