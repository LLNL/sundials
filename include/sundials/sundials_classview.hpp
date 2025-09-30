/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * Base classes for C++ implementations wrappers (views) of SUNDIALS objects.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CLASSVIEW_HPP
#define _SUNDIALS_CLASSVIEW_HPP

#include <sundials/sundials_convertibleto.hpp>
#include <utility>

namespace sundials {
namespace experimental {

template<class T, class Deleter>
class ClassView : public sundials::ConvertibleTo<T>
{
public:
  ClassView() : object_(nullptr) {}

  ClassView(T& object) : object_(object) {}

  ClassView(T&& object) : object_(std::forward<T>(object)) {}

  ClassView(const ClassView&) = delete;

  ClassView(ClassView&& other) noexcept
    : object_(std::exchange(other.object_, nullptr))
  {}

  ClassView& operator=(const ClassView&) = delete;

  ClassView& operator=(ClassView&& rhs) noexcept
  {
    this->object_ = std::exchange(rhs.object_, nullptr);
    return *this;
  };

  ~ClassView()
  {
    if (object_) { Deleter{}(this->get()); }
  };

  // Override ConvertibleTo functions
  T get() override { return object_; }

  T get() const override { return object_; }

  operator T() override { return object_; }

  operator T() const override { return object_; }

protected:
  T object_;
};

template<typename T, typename Deleter, typename Func, typename... Args>
T Create(Args&&... args)
{
  return ClassView<T, Deleter>(Func(std::forward<Args>(args)...));
}

} // namespace experimental
} // namespace sundials

#endif
