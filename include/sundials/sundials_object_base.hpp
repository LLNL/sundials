/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
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
 *
 * Generic support for C++ classes that wrap SUNDIALS C objects.
 *
 * INTERNAL USE: this header is the generic layer under the type-specific bases
 * (e.g., sundials::NVectorBase in sundials_nvector_base.hpp). It is installed
 * because those headers include it, but its contents are in namespace
 * sundials::impl and are not a user interface. It provides
 *
 * - error reporting at the C boundary: handleActiveException and the guards
 *   guardVoid, guardErr, guardBool and guardValue;
 * - override detection: impl::is_overridden_v and impl::same_signature_v;
 * - the lifetime layer impl::ObjectBase<Derived, Traits, TypeBase>.
 *
 * Error reporting notes
 *
 * - Errors are reported with the file and line passed by the caller (each
 *   trampoline passes __FILE__ and __LINE__), so the error handler's output
 *   points at the trampoline rather than at the shared helper.
 * - handleActiveException and the guards report through the public C function
 *   SUNContext_HandleError, which logs the message verbatim: a % in the message
 *   is not treated as a format directive and the message is not truncated.
 * - A NULL message is replaced by the default text of the error code. A NULL
 *   SUNContext goes to the global fallback handler, which prints to stderr.
 * - An error handler that throws, when called from a noexcept function such as
 *   a trampoline, ends in std::terminate.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_SUNDIALS_OBJECT_BASE_HPP
#define SUNDIALS_SUNDIALS_OBJECT_BASE_HPP

#if defined(_MSVC_LANG)
#if _MSVC_LANG < 201703L
#error "sundials_object_base.hpp requires C++17 or later"
#endif
#elif !defined(__cplusplus) || __cplusplus < 201703L
#error "sundials_object_base.hpp requires C++17 or later"
#endif

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

namespace sundials {
namespace impl {

// -----------------------------------------------------------------------------
// Error reporting at the C boundary
// -----------------------------------------------------------------------------

// Classifies the exception currently being handled, reports it, and returns the
// error code. Call it only from inside a catch block: SUNContext_HandleError is
// called while the matching handler is active, so what() is still valid.
//
// - std::bad_alloc gives SUN_ERR_MEM_FAIL
// - any other std::exception gives SUN_ERR_OP_FAIL with what() as the message
// - anything else gives SUN_ERR_OP_FAIL
inline SUNErrCode handleActiveException(SUNContext ctx, const char* func,
                                        const char* file, int line) noexcept
{
  try
  {
    throw;
  }
  catch (const std::bad_alloc& err)
  {
    SUNContext_HandleError(ctx, SUN_ERR_MEM_FAIL, line, func, file, err.what());
    return SUN_ERR_MEM_FAIL;
  }
  catch (const std::exception& err)
  {
    SUNContext_HandleError(ctx, SUN_ERR_OP_FAIL, line, func, file, err.what());
    return SUN_ERR_OP_FAIL;
  }
  catch (...)
  {
    SUNContext_HandleError(ctx, SUN_ERR_OP_FAIL, line, func, file,
                           "unknown C++ exception");
    return SUN_ERR_OP_FAIL;
  }
}

// Guards: run callable and catch every exception so it cannot propagate into C code.
// They are named by the kind of operation because SUNErrCode and sunbooleantype
// are both int. Each takes the C API name and the file and line of the caller
// (__FILE__, __LINE__), which are reported with the error.

// For operations returning void
template<class Callable>
void guardVoid(SUNContext ctx, const char* func, const char* file, int line,
               Callable&& callable) noexcept
{
  try
  {
    callable();
  }
  catch (...)
  {
    handleActiveException(ctx, func, file, line);
  }
}

// For operations returning SUNErrCode: returns the error code on an exception
template<class Callable>
SUNErrCode guardErr(SUNContext ctx, const char* func, const char* file,
                    int line, Callable&& callable) noexcept
{
  try
  {
    return callable();
  }
  catch (...)
  {
    return handleActiveException(ctx, func, file, line);
  }
}

// For operations returning sunbooleantype: returns SUNFALSE on an exception
template<class Callable>
sunbooleantype guardBool(SUNContext ctx, const char* func, const char* file,
                         int line, Callable&& callable) noexcept
{
  try
  {
    return callable();
  }
  catch (...)
  {
    handleActiveException(ctx, func, file, line);
  }
  return SUNFALSE;
}

// For all other operations: returns placeholder on an exception only because
// the C signature requires a value. It is not an error indicator.
template<class Result, class Callable>
Result guardValue(SUNContext ctx, const char* func, const char* file, int line,
             Result placeholder, Callable&& callable) noexcept
{
  try
  {
    return callable();
  }
  catch (...)
  {
    handleActiveException(ctx, func, file, line);
  }
  return placeholder;
}

// -----------------------------------------------------------------------------
// Override detection
//
// For a member function name X, &Derived::X names the X found by lookup starting
// in Derived. If Derived does not declare X, the expression has the type of the
// base class declaration, e.g., void (Base::*)(const Derived&). If Derived (or
// an intermediate class) declares X, the class part of the type differs. The
// traits below compare the two types:
//
//   is_overridden_v<DerivedFn, BaseFn>   true if the two types differ, i.e.,
//                                        Derived redeclares the function
//   same_signature_v<DerivedFn, BaseFn>  true if both are member functions with
//                                        the same return type, parameters and
//                                        constness (the class part is ignored)
//
// In C++17 noexcept is part of the function type, so there is one
// specialization for each combination of const and noexcept, and an override
// that adds noexcept matches the base declaration.
// -----------------------------------------------------------------------------

// Splits a pointer to member function type. For anything else (a data member
// or a static function with the same name) is_member_fn is false, so the caller
// can report a readable error instead of a hard compile error.
template<class MemberFn>
struct member_fn_signature
{
  static constexpr bool is_member_fn = false;
};

template<class Result, class Class, class... Args>
struct member_fn_signature<Result (Class::*)(Args...)>
{
  static constexpr bool is_member_fn = true;
  static constexpr bool is_const     = false;
  using type                         = Result(Args...);
};

template<class Result, class Class, class... Args>
struct member_fn_signature<Result (Class::*)(Args...) const>
{
  static constexpr bool is_member_fn = true;
  static constexpr bool is_const     = true;
  using type                         = Result(Args...);
};

template<class Result, class Class, class... Args>
struct member_fn_signature<Result (Class::*)(Args...) noexcept>
{
  static constexpr bool is_member_fn = true;
  static constexpr bool is_const     = false;
  using type                         = Result(Args...);
};

template<class Result, class Class, class... Args>
struct member_fn_signature<Result (Class::*)(Args...) const noexcept>
{
  static constexpr bool is_member_fn = true;
  static constexpr bool is_const     = true;
  using type                         = Result(Args...);
};

template<class DerivedFn, class BaseFn>
inline constexpr bool is_overridden_v = !std::is_same_v<DerivedFn, BaseFn>;

template<class DerivedFn, class BaseFn,
         bool = member_fn_signature<DerivedFn>::is_member_fn&&
           member_fn_signature<BaseFn>::is_member_fn>
struct same_signature : std::false_type
{};

template<class DerivedFn, class BaseFn>
struct same_signature<DerivedFn, BaseFn, true>
  : std::bool_constant<
      std::is_same_v<typename member_fn_signature<DerivedFn>::type,
                     typename member_fn_signature<BaseFn>::type> &&
      member_fn_signature<DerivedFn>::is_const ==
        member_fn_signature<BaseFn>::is_const>
{};

template<class DerivedFn, class BaseFn>
inline constexpr bool same_signature_v = same_signature<DerivedFn, BaseFn>::value;

// -----------------------------------------------------------------------------
// ObjectBase
//
// Generic lifetime base for C++ wrappers of SUNDIALS C objects that follow the
// "NewEmpty" pattern: a struct with `content` and `sunctx` fields that is
// created by a NewEmpty function and freed by a FreeEmpty function (which free
// the struct and its ops table but not the content).
//
// ObjectBase owns the handle (obj_handle_) and the obj_c_owned_ flag, and provides
// the special members, the handle accessors, and the clone and destroy helpers.
// It only touches `content` and `sunctx`: filling the ops table belongs to the
// type-specific base (the TypeBase template argument, e.g., NVectorBase<Derived>).
//
// Hook contract. TypeBase provides a private static
//
//   static void sunInstallOps(Struct*) noexcept
//
// that ObjectBase calls on every new handle: the struct is allocated and zeroed
// by NewEmpty, content points at the ObjectBase subobject, and the hook only
// fills the ops table. TypeBase declares ObjectBase as a friend so that the hook
// can stay private. The call is qualified (TypeBase::sunInstallOps), so a
// member of Derived with the same name neither hides nor replaces the
// installer. Derived must inherit publicly from TypeBase and pass itself as the
// template argument; both are checked, with one error message.
//
// Handle storage. content holds the address of the ObjectBase subobject, not
// of Derived: the handle is created in a base class constructor, before Derived
// exists, so a downcast there would be undefined behavior. objFromHandle() does the
// downcast in two steps, void* to ObjectBase* to Derived*, which applies the
// correct offset when ObjectBase is not at the start of Derived (a polymorphic
// first base, or multiple inheritance).
//
// Ownership. The flag obj_c_owned_ is set only by objAdoptClone, which a clone
// operation calls on an object it allocated with new. objDestroy, the body
// of the C destroy op, deletes an object only if obj_c_owned_ is set; for a
// user-owned object it does nothing and the object's destructor frees the
// handle. The destroy op must always be installed, so the C fallback that
// frees content never runs. Each type-specific base supplies a trampoline with
// its type's destroy signature that calls objDestroy.
//
// Special members (all private, used by the type-specific base, whose own
// special members a derived class normally defaults):
//
//   operation         source        target before  result
//   default ctor      -             -              no handle
//   copy ctor         no handle     -              no handle
//   copy ctor         has handle    -              new handle from the source's
//                                                  context; the source is
//                                                  unchanged (may throw)
//   move ctor         no handle     -              no handle
//   move ctor         non-clone     -              takes the handle and repoints
//                     with handle                  content; source has no handle
//   move ctor         clone         -              new handle from the source's
//                                                  context; the clone keeps its
//                                                  handle; on failure the error
//                                                  is reported and the new
//                                                  object has no handle
//   copy assign       any           has handle     keeps its handle
//   copy assign       has handle    no handle      new handle from the source's
//                                                  context (may throw)
//   copy assign       no handle     no handle      no handle
//   move assign       any           has handle     keeps its handle; the source
//                                                  keeps its own
//   move assign       any           no handle      as the move constructor
//   self-assignment   -             -              no change
//
// obj_c_owned_ is never copied, moved or transferred. The default constructor and
// both moves are noexcept and no path terminates; the copy members may throw.
//
// Template parameters:
//   Derived  the concrete C++ type (CRTP)
//   TypeBase the type-specific base class that provides the hook
//   Traits   per-type details, which must provide
//              using Struct = ...;  a struct with `content` and `sunctx`
//              static Struct* NewEmpty(SUNContext);  NULL on failure
//              static void FreeEmpty(Struct*);
//              static constexpr const char* CLASS_NAME = "...";  for messages
//
// Why a traits type instead of NewEmpty and FreeEmpty as template arguments:
//
// - Windows shared builds. A function pointer template argument must be a
//   constant expression. With SUNDIALS_EXPORT as dllimport, the address of
//   N_VNewEmpty is not one on GCC/MinGW (it is resolved through the import
//   table at load time), so ObjectBase<..., &N_VNewEmpty, ...> fails to compile
//   there. The inline wrappers in a traits type are ordinary header functions
//   that call the imported function, so they work everywhere.
// - Readability. ObjectBase<MyVec, NVectorTraits, NVectorBase<MyVec>> replaces a list of function
//   pointers, and the traits are the one place for per-type details (struct
//   type, create and free functions, the class name used in messages).
// -----------------------------------------------------------------------------

template<class Derived, class Traits, class TypeBase>
class ObjectBase
{
  // All members are private. Only the type-specific base (TypeBase, e.g.,
  // NVectorBase<Derived>) may use them: Derived may not, because calling
  // objRelease or objAdoptClone from a user class breaks ownership, and the obj
  // prefix keeps the names out of the way of members of Derived. TypeBase in
  // turn befriends ObjectBase so that ObjectBase can call its private
  // sunInstallOps. Both friend declarations are needed.
  friend TypeBase;

private:
  using ObjStruct = typename Traits::Struct;

  // ---------------------------------------------------------------------------
  // Construction
  //
  // The default constructor creates an object without a handle. The SUNContext
  // constructor allocates the handle with NewEmpty (which nullifies the content
  // and all ops pointers), sets the content pointer to this object, and
  // installs the ops through the per-type hook. It throws std::invalid_argument
  // if ctx is NULL and std::bad_alloc if the handle cannot be allocated.
  // ---------------------------------------------------------------------------

  ObjectBase() noexcept = default;

  explicit ObjectBase(SUNContext ctx) { objAllocate(ctx); }

  // ---------------------------------------------------------------------------
  // Copy and move
  //
  // - Copy construction: an empty source gives an empty object. Otherwise a new
  //   handle is allocated from the source's context (may throw). The source is
  //   unchanged.
  // - Move construction: an empty source gives an empty object. A source that
  //   is not a clone gives up its handle (content is repointed to the new
  //   object) and becomes empty. A clone (obj_c_owned_) keeps its handle, because
  //   SUNDIALS holds it, so a new handle is allocated from the source's
  //   context. If that fails, the error is reported through the context and the
  //   new object is empty (the derived class still moves its data).
  // - Assignment never removes the target's handle: a target with a handle
  //   keeps it and only the derived class's data is assigned (also from an
  //   empty source). An empty target acts like the corresponding constructor.
  // - obj_c_owned_ is never copied, moved or transferred. It is set only by
  //   objAdoptClone.
  // - Moves and the default constructor are noexcept; the copy members may
  //   throw.
  // ---------------------------------------------------------------------------

  ObjectBase(const ObjectBase& other)
  {
    if (other.obj_handle_) { objAllocate(other.objContext()); }
  }

  ObjectBase(ObjectBase&& other) noexcept { objMoveHandleFrom(other); }

  ObjectBase& operator=(const ObjectBase& other)
  {
    if (this != &other && obj_handle_ == nullptr && other.obj_handle_ != nullptr) { objAllocate(other.objContext()); }
    return *this;
  }

  ObjectBase& operator=(ObjectBase&& other) noexcept
  {
    if (this != &other && obj_handle_ == nullptr) { objMoveHandleFrom(other); }
    return *this;
  }

  // The destructor is private and non-virtual: deleting through an
  // ObjectBase* does not compile. Its body is exactly objRelease(). Clones reach it
  // through the delete in objDestroy; user-owned objects through their
  // normal destruction.
  ~ObjectBase() { objRelease(); }

  // The handle (NULL if there is none).
  ObjStruct* objHandle() const noexcept { return obj_handle_; }

  bool objHasHandle() const noexcept { return obj_handle_ != nullptr; }

  // The context of the handle (NULL if there is no handle). Reading it from the
  // handle keeps it consistent with SUNDIALS, which may overwrite it.
  SUNContext objContext() const noexcept { return obj_handle_ ? obj_handle_->sunctx : nullptr; }

  // The object for a handle. The handle stores the ObjectBase subobject, so
  // this is a two-step cast: with multiple inheritance or a polymorphic first
  // base, the ObjectBase subobject is not at the address of Derived.
  static Derived* objFromHandle(ObjStruct* handle) noexcept
  {
    return static_cast<Derived*>(static_cast<ObjectBase*>(handle->content));
  }

  // Allocates the handle (obj_handle_ must be NULL). Throws as the constructor does.
  void objAllocate(SUNContext ctx)
  {
    if (ctx == nullptr)
    {
      throw std::invalid_argument(std::string(Traits::CLASS_NAME) +
                                  ": SUNContext is NULL");
    }
    ObjStruct* handle = objNewHandle(ctx);
    if (handle == nullptr) { throw std::bad_alloc(); }
    obj_handle_ = handle;
  }

  // Detaches the content pointer and frees the handle. Only the content is
  // detached: it belongs to the C++ object, and FreeEmpty frees the ops table
  // and the struct itself. After release the object has no handle.
  void objRelease() noexcept
  {
    if (obj_handle_)
    {
      obj_handle_->content = nullptr;
      Traits::FreeEmpty(obj_handle_);
      obj_handle_ = nullptr;
    }
  }

  // Marks an object allocated with new by a clone operation as owned by its
  // handle (its lifetime is then controlled by the C destroy op) and returns
  // the handle.
  static ObjStruct* objAdoptClone(Derived* clone) noexcept
  {
    ObjectBase* self = clone;
    self->obj_c_owned_ = true;
    return self->obj_handle_;
  }

  // The body of the C destroy op (always installed). Only objects marked by
  // objAdoptClone are deleted here; the destructor then frees the handle, so there
  // is a single free path. For any other object, which the user owns, this does
  // nothing and the object's destructor frees the handle.
  //
  // It deletes through Derived*, so Derived's destructor must be public (or
  // Derived must befriend this class).
  static void objDestroy(ObjStruct* handle) noexcept
  {
#if defined(SUNDIALS_DEBUG)
    if (handle->content == nullptr)
    {
      std::fprintf(stderr, "%s:%d: %s: ObjectBase: destroy called on a handle "
                           "with NULL content\n",
                   __FILE__, __LINE__, __func__);
      std::abort();
    }
#endif
    ObjectBase* self = static_cast<ObjectBase*>(handle->content);
    if (self->obj_c_owned_) { delete static_cast<Derived*>(self); }
  }

  // Creates a handle for this object and installs the ops. Returns NULL on
  // failure and reports nothing.
  ObjStruct* objNewHandle(SUNContext ctx) noexcept
  {
    ObjStruct* handle = Traits::NewEmpty(ctx);
    if (handle == nullptr) { return nullptr; }

    // content holds the ObjectBase subobject, not Derived: Derived is not
    // constructed yet when this runs from a constructor.
    handle->content = static_cast<void*>(this);

    // True if Derived derives from this class and the conversion is accessible
    // (public inheritance). Both the downcast in objFromHandle and the trampolines
    // need it. Local, so it is evaluated only when Derived is complete.
    constexpr bool derives_publicly =
      std::is_base_of_v<ObjectBase, Derived> &&
      std::is_convertible_v<Derived*, ObjectBase*>;

    static_assert(derives_publicly,
                  "sundials: Derived must inherit publicly from the SUNDIALS "
                  "C++ base class (e.g., sundials::NVectorBase<Derived>) and "
                  "pass itself as the template argument");

    // Qualified call: it cannot be replaced by a member of Derived.
    TypeBase::sunInstallOps(handle);
    return handle;
  }

  // Takes the handle of other, or allocates one if other is a clone. Only for an object
  // without a handle.
  void objMoveHandleFrom(ObjectBase& other) noexcept
  {
    if (other.obj_handle_ == nullptr) { return; }
    if (!other.obj_c_owned_)
    {
      obj_handle_          = other.obj_handle_;
      other.obj_handle_        = nullptr;
      obj_handle_->content = static_cast<void*>(this);
    }
    else
    {
      SUNContext ctx = other.objContext();
      obj_handle_             = objNewHandle(ctx);
      if (obj_handle_ == nullptr)
      {
        SUNContext_HandleError(ctx, SUN_ERR_MEM_FAIL, __LINE__,
                               Traits::CLASS_NAME, __FILE__,
                               "out of memory allocating a handle while moving "
                               "from a clone; the new object has no handle");
      }
    }
  }

  // The C struct (NULL means the object has no handle).
  ObjStruct* obj_handle_ = nullptr;

  // True if this object was allocated by a clone operation and its lifetime is
  // controlled by the C destroy op. Set only by objAdoptClone; it describes this
  // object's storage, not the handle, and is never copied, moved or transferred.
  bool obj_c_owned_ = false;
};

} // namespace impl
} // namespace sundials

#endif // SUNDIALS_SUNDIALS_OBJECT_BASE_HPP
