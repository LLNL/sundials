/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * C++ interface to the SUNDIALS context object
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_HPP
#define _SUNDIALS_CONTEXT_HPP

#include <memory>
#include <sundials/sundials_context.h>
#include <sundials/sundials_convertibleto.hpp>
#include "sundials/sundials_types.h"

namespace sundials {

class Context : public sundials::ConvertibleTo<SUNContext>
{
public:
  explicit Context(SUNComm comm = SUN_COMM_NULL)
  {
    sunctx_ = std::make_unique<SUNContext>();
    SUNContext_Create(comm, sunctx_.get());
  }

  /* disallow copy, but allow move construction */
  Context(const Context&) = delete;
  Context(Context&&)      = default;

  /* disallow copy, but allow move operators */
  Context& operator=(const Context&) = delete;
  Context& operator=(Context&&) = default;

  SUNContext Convert() override
  {
    return *sunctx_.get();
  }
  SUNContext Convert() const override
  {
    return *sunctx_.get();
  }
  operator SUNContext() override
  {
    return *sunctx_.get();
  }
  operator SUNContext() const override
  {
    return *sunctx_.get();
  }

  ~Context()
  {
    if (sunctx_) SUNContext_Free(sunctx_.get());
  }

private:
  std::unique_ptr<SUNContext> sunctx_;
};

} // namespace sundials

#endif // _SUNDIALS_CONTEXT_HPP
