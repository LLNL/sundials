 /* -----------------------------------------------------------------
  * SUNDIALS Copyright Start
  * Copyright (c) 2002-2022, Lawrence Livermore National Security
  * and Southern Methodist University.
  * All rights reserved.
  *
  * See the top-level LICENSE and NOTICE files for details.
  *
  * SPDX-License-Identifier: BSD-3-Clause
  * SUNDIALS Copyright End
  * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_HPP
#define _SUNDIALS_CONTEXT_HPP

#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_context.h>

namespace sundials {

// RAII wrapper of SUNContext.
// NOTE: SUNContext_Free should not be called on the underlying or a double free will occur!
class Context : public impl::ClassView<_SUNContext, std::unique_ptr<_SUNContext>> {
public:
  explicit Context(void* comm = nullptr);
  virtual ~Context() override;
private:
  SUNContext sunctx;
};

} // namespace sundials

#endif // _SUNDIALS_CONTEXT_HPP
