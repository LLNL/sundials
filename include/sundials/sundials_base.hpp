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
 * Base classes for C++ implementations of SUNDIALS objects and wrappers (views)
 * of SUNDIALS objects
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_BASE_HPP
#define _SUNDIALS_BASE_HPP

#include <memory>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_convertibleto.hpp>

namespace sundials {
namespace impl {

//
// Common base class for C++ implementations of SUNDIALS data structures.
//
template<class ObjectStruct, class ObjectOps>
class BaseObject
{
public:
  BaseObject() = default;

  BaseObject(SUNContext sunctx)
    : sunctx_(sunctx),
      object_(std::make_unique<ObjectStruct>()),
      object_ops_(std::make_unique<ObjectOps>())
  {
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
  }

  // Move constructor
  BaseObject(BaseObject&& other) noexcept
    : sunctx_(std::move(other.sunctx_)),
      object_(std::move(other.object_)),
      object_ops_(std::move(other.object_ops_))
  {
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
  }

  // Copy constructor
  BaseObject(const BaseObject& other)
    : sunctx_(other.sunctx_),
      object_(std::make_unique<ObjectStruct>()),
      object_ops_(std::make_unique<ObjectOps>(*other.object_ops_))
  {
    object_->content = this;
    object_->sunctx  = other.sunctx_;
    object_->ops     = object_ops_.get();
  }

  // Move assignment
  BaseObject& operator=(BaseObject&& rhs) noexcept
  {
    sunctx_          = std::move(rhs.sunctx_);
    object_ops_      = std::move(rhs.object_ops_);
    object_          = std::move(rhs.object_);
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
    return *this;
  }

  // Copy assignment
  BaseObject& operator=(const BaseObject& rhs)
  {
    sunctx_          = rhs.sunctx_;
    object_ops_      = std::make_unique<ObjectOps>(*rhs.object_ops_);
    object_          = std::make_unique<ObjectStruct>();
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
    return *this;
  }

  // We have a pure virtual destructor to make this an asbtract class
  virtual ~BaseObject() = 0;

  // Getters
  SUNContext sunctx() const { return this->object_->sunctx; }

protected:
  // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
  SUNContext sunctx_{};
  // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
  std::unique_ptr<ObjectStruct> object_;
  // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
  std::unique_ptr<ObjectOps> object_ops_;
};

// Pure virtual destructor requires implementation
template<class ObjectStruct, class ObjectOps>
BaseObject<ObjectStruct, ObjectOps>::~BaseObject() = default;

} // namespace impl

namespace experimental {

template<class T, class Deleter>
class ClassView : public sundials::ConvertibleTo<T>
{
public:
  ClassView() : object_(nullptr) {}

  ClassView(T&& object) : object_(std::make_unique<T>(object)) {}

  ClassView(const ClassView&)  = delete;
  ClassView(ClassView&& other) = default;

  ClassView& operator=(const ClassView&) = delete;
  ClassView& operator=(ClassView&& rhs)  = default;

  ~ClassView()
  {
    if (object_) { Deleter{}(this->Convert()); }
  };

  // Override ConvertibleTo functions
  T Convert() override { return *object_.get(); }

  T Convert() const override { return *object_.get(); }

  operator T() override { return *object_.get(); }

  operator T() const override { return *object_.get(); }

private:
  std::unique_ptr<T> object_;
};

} // namespace experimental
} // namespace sundials

#endif
