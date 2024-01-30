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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_PROFILER_HPP
#define _SUNDIALS_PROFILER_HPP

#include <cstring>
#include <sundials/sundials_config.h>
#include <sundials/sundials_profiler.h>

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && defined(SUNDIALS_CALIPER_ENABLED)
#define SUNDIALS_CXX_MARK_FUNCTION(projobj) CALI_CXX_MARK_FUNCTION
#elif defined(SUNDIALS_BUILD_WITH_PROFILING)
#define SUNDIALS_CXX_MARK_FUNCTION(profobj) \
  sundials::ProfilerMarkScope ProfilerMarkScope__(profobj, __func__)
#else
#define SUNDIALS_CXX_MARK_FUNCTION(profobj)
#endif

namespace sundials {
/* Convenience class for C++ codes.
   Allows for simpler profiler statements using C++ scoping rules. */
class ProfilerMarkScope
{
public:
  ProfilerMarkScope(SUNProfiler prof, const char* name)
  {
    prof_ = prof;
    name_ = name;
    SUNProfiler_Begin(prof_, name_);
  }

  ~ProfilerMarkScope() { SUNProfiler_End(prof_, name_); }

private:
  SUNProfiler prof_;
  const char* name_;
};
} // namespace sundials

#endif /* SUNDIALS_PROFILER_HPP */
