/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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

#ifndef _SUNDIALS_LOGGER_HPP
#define _SUNDIALS_LOGGER_HPP

#include <sundials/sundials_config.h>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_logger.h>

namespace sundials {
namespace experimental {

class SUNLoggerView : public sundials::ConvertibleTo<SUNLogger>
{
public:
  SUNLoggerView(SUNComm comm, int output_rank)
  {
    logger_ = std::make_unique<SUNLogger>();
    SUNLogger_Create(comm, output_rank, logger_.get());
  }

  SUNLoggerView(SUNComm comm)
  {
    logger_ = std::make_unique<SUNLogger>();
    SUNLogger_CreateFromEnv(comm, logger_.get());
  }

  SUNLoggerView(SUNLogger sunlogger) 
  {  logger_.reset(&sunlogger);  }

  /* disallow copy, but allow move construction */
  SUNLoggerView(const SUNLoggerView&) = delete;
  SUNLoggerView(SUNLoggerView&&)      = default;

  /* disallow copy, but allow move operators */
  SUNLoggerView& operator=(const SUNLoggerView&) = delete;
  SUNLoggerView& operator=(SUNLoggerView&&)      = default;

  SUNLogger get() override { return *logger_.get(); }

  SUNLogger get() const override { return *logger_.get(); }

  operator SUNLogger() override { return *logger_.get(); }

  operator SUNLogger() const override { return *logger_.get(); }

  template<typename... Args>
  static SUNLoggerView Create(Args&&... args);

  ~SUNLoggerView()
  {
    if (logger_) { SUNLogger_Destroy(logger_.get()); }
  }

private:
  std::unique_ptr<SUNLogger> logger_;
};


template<typename... Args>
SUNLoggerView SUNLoggerView::Create(Args&&... args)
{
  return SUNLoggerView(std::forward<Args>(args)...);
}


} // namespace experimental
} // namespace sundials

#endif /* SUNDIALS_LOGGER_HPP */
