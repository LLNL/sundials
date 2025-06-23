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

#ifndef _PYSUNDIALS_ARKODE_USERSUPPLIED_HPP
#define _PYSUNDIALS_ARKODE_USERSUPPLIED_HPP

#include "pysundials_helpers.hpp"

////////////////////////////////////////////////////////////////////////////////
// std::function types corresponding to C function pointer typedefs
///////////////////////////////////////////////////////////////////////////////

using ARKRhsStdFn = int(sunrealtype, N_Vector, N_Vector, void*);

///////////////////////////////////////////////////////////////////////////////
// ERKStep user-supplied functions
///////////////////////////////////////////////////////////////////////////////

struct erkstep_user_supplied_fn_table
{
  nb::object erkstep_rhsfn;
};

inline int erk_rhsfn_wrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                             void* user_data)
{
  return pysundials::user_supplied_fn_wrapper<
    ARKRhsStdFn,
    erkstep_user_supplied_fn_table>(t, y, ydot, user_data,
                                    &erkstep_user_supplied_fn_table::erkstep_rhsfn);
}

#endif
