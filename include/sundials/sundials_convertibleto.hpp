/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * Base class for converting C++ wrappers (views) to SUNDIALS objects
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CONVERTIBLETO_HPP
#define _SUNDIALS_CONVERTIBLETO_HPP

namespace sundials {

template<class T>
class ConvertibleTo
{
public:
  // Explicit conversion to the underlying type
  virtual T Convert()       = 0;
  virtual T Convert() const = 0;

  // Implicit conversion to the underlying type
  virtual operator T()       = 0;
  virtual operator T() const = 0;

  virtual ~ConvertibleTo() = default;
};

} // namespace sundials

#endif // _SUNDIALS_CONVERTIBLETO_HPP
