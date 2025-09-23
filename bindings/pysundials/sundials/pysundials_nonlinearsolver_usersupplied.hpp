/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef _PYSUNDIALS_NONLINEARSOLVER_USERSUPPLIED_HPP
#define _PYSUNDIALS_NONLINEARSOLVER_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nonlinearsolver.hpp>

#include "pysundials_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

struct SUNNonlinearSolverFunctionTable
{
  nb::object sysfn;
  nb::object lsetupfn;
  nb::object lsolvefn;
  nb::object convtestfn;
};

inline SUNNonlinearSolverFunctionTable* SUNNonlinearSolverFunctionTable_Alloc()
{
  auto fn_table = static_cast<SUNNonlinearSolverFunctionTable*>(
    std::malloc(sizeof(SUNNonlinearSolverFunctionTable)));
  std::memset(fn_table, 0, sizeof(SUNNonlinearSolverFunctionTable));
  return fn_table;
}

template<typename... Args>
inline int sunnonlinearsolver_sysfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolSysFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::sysfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_lsetupfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolLSetupFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::lsetupfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_lsolvefn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolLSolveFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::lsolvefn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_convtestfn_wrapper(Args... args)
{
  return pysundials::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolConvTestFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::convtestfn, std::forward<Args>(args)...);
}

#endif // _PYSUNDIALS_NONLINEARSOLVER_USERSUPPLIED_HPP
