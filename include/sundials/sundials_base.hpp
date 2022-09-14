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

#ifndef _SUNDIALS_BASE_HPP
#define _SUNDIALS_BASE_HPP

#include <memory>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>

namespace sundials {

template<class T>
class ConvertibleTo {
public:
  // Explicit conversion to the underlying type
  virtual T get() = 0;
  virtual T get() const = 0;

  // Implicit conversion to the underlying type
  virtual operator T() = 0;
  virtual operator T() const = 0;

  virtual ~ConvertibleTo() = default;
};

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
      : sunctx_(sunctx), object_(std::make_unique<ObjectStruct>()), object_ops_(std::make_unique<ObjectOps>())
  {
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
  }

  // Move constructor
  BaseObject(BaseObject&& other) noexcept
      : sunctx_(std::move(other.sunctx_)), object_(std::move(other.object_)), object_ops_(std::move(other.object_ops_))
  {
    object_->content = this;
    object_->sunctx  = sunctx_;
    object_->ops     = object_ops_.get();
  }

  // Copy constructor clones the gko::matrix and SUNMatrix
  BaseObject(const BaseObject& other)
      : sunctx_(other.sunctx_), object_(std::make_unique<ObjectStruct>()),
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

  // Copy assignment clones the gko::matrix and SUNMatrix
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

template<class T, class PtrType = std::unique_ptr<T>, class Deleter = std::default_delete<T*>>
class ClassView : public ConvertibleTo<T*> {
public:
  ClassView<T, PtrType, Deleter>()
    : underlying_ptr_(new T()) {};

  // Only allow construction from Rvalue since we take complete ownership
  ClassView<T, PtrType, Deleter>(T*&& val) : underlying_ptr_(val) {}
  ClassView<T, PtrType, Deleter>(T*& val) = delete;

  // Disallow copy constructor, use default move constructor
  ClassView<T, PtrType, Deleter>(const ClassView<T, PtrType, Deleter>&) = delete;
  ClassView<T>(ClassView<T, PtrType, Deleter>&&) noexcept = default;

  // Disallow copy assignment, use default move assignment
  ClassView<T, PtrType, Deleter>& operator=(const ClassView<T, PtrType, Deleter>&) = delete;
  ClassView<T, PtrType, Deleter>& operator=(ClassView<T, PtrType, Deleter>&&) noexcept = default;

  // ClassView is an abstract class
  virtual ~ClassView<T, PtrType, Deleter>() = 0;

  // Override ConvertibleTo functions
  virtual T* get() override { return underlying_ptr_.get(); }
  virtual T* get() const override { return underlying_ptr_.get(); }
  virtual operator T*() override { return underlying_ptr_.get(); }
  virtual operator T*() const override { return underlying_ptr_.get(); }

protected:
  PtrType underlying_ptr_;
};

// Pure virtual destructor requires a definition.
template<class T, class PtrType, class Deleter>
ClassView<T, PtrType, Deleter>::~ClassView<T, PtrType, Deleter> () = default;

} // namespace impl
} // namespace sundials

#endif
