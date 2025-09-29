/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_HPP
#define _SUNDIALS_LOGGER_HPP

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_logger.h>

namespace sundials {
namespace experimental {

struct SUNLoggerDeleter
{
  void operator()(SUNLogger logger) { SUNLogger_Destroy(&logger); }
};

class SUNLoggerView
  : public sundials::experimental::ClassView<SUNLogger, SUNLoggerDeleter>
{
public:
  using sundials::experimental::ClassView<SUNLogger, SUNLoggerDeleter>::ClassView;

  SUNLoggerView(SUNComm comm, int output_rank)
  {
    SUNLogger_Create(comm, output_rank, &object_);
  }

  SUNLoggerView(SUNComm comm) { SUNLogger_CreateFromEnv(comm, &object_); }

  template<typename... Args>
  static SUNLoggerView Create(Args&&... args);
};

template<typename... Args>
SUNLoggerView SUNLoggerView::Create(Args&&... args)
{
  return SUNLoggerView(std::forward<Args>(args)...);
}

} // namespace experimental
} // namespace sundials

#endif /* SUNDIALS_LOGGER_HPP */
