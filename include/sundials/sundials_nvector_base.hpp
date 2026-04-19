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
 * This header provides the NVectorBase<Derived> base class that bridges C++
 * vector implementations to the SUNDIALS C N_Vector interface.
 *
 * -----------------------------------------------------------------------------
 * OVERVIEW
 * -----------------------------------------------------------------------------
 *
 * The three key aspects of the bridge class design are:
 *
 * 1. The derived class inherits from NVectorBase<Derived>, giving the base
 *    class compile-time knowledge of the concrete type, so operations can take
 *    and return the derived type (e.g., const MyVec&).
 *
 * 2. Every N_Vector operation is a virtual member function of the base class.
 *    The derived class overrides the operations it implements. A small core of
 *    operations is required (pure virtual), and an operation that the derived
 *    class does not override is not installed in the ops table (see OPERATIONS
 *    below).
 *
 * 3. For each installed operation, the base class defines a static function
 *    with the C signature SUNDIALS expects ("trampoline"). Each trampoline
 *    extracts the C++ objects from the N_Vectors and calls the virtual
 *    function.
 *
 * Derived class methods are named after the N_Vector operation they implement
 * without the N_V prefix, e.g., LinearSum implements N_VLinearSum.
 *
 * -----------------------------------------------------------------------------
 * EXAMPLE
 * -----------------------------------------------------------------------------
 *
 * class MyVec final : public sundials::NVectorBase<MyVec> // final recommended
 * {
 * public:
 *   MyVec() = default;
 *
 *   MyVec(sunindextype length, SUNContext ctx)
 *     : NVectorBase(ctx), vec_(static_cast<std::size_t>(length), 0.0)
 *   {}
 *
 *   // Used by N_VClone. The new object must be allocated with new and its
 *   // values do not need to be initialized.
 *   MyVec* Clone(SUNContext ctx) const override
 *   {
 *     return new MyVec(GetLength(), ctx);
 *   }
 *
 *   sunindextype GetLength() const override
 *   {
 *     return static_cast<sunindextype>(vec_.size());
 *   }
 *
 *   // z = a*x + b*y, where z is *this. x and/or y may be the same object as
 *   // *this (see ALIASING below).
 *   void LinearSum(sunrealtype a, const MyVec& x, sunrealtype b,
 *                  const MyVec& y) override;
 *
 *   // Const, Scale, DotProd, MaxNorm (required)
 *   // ...
 *
 *   // Optional operations that the packages and features in use need (see
 *   // OPERATIONS below)
 *   void Div(const MyVec& x, const MyVec& y) override;
 *   sunrealtype WrmsNorm(const MyVec& w) const override;
 *   SUNErrCode LinearCombination(ConstReals c, ConstVectors X) override;
 *
 * private:
 *   std::vector<sunrealtype> vec_;
 * };
 *
 * Rules for the derived class:
 *
 * - Inherit publicly from NVectorBase<Derived> and pass the derived class
 *   itself as the template argument. Both are checked at compile time.
 *
 * - Use override on every operation. Without it, a method whose signature does
 *   not exactly match the base class declaration becomes a new function. For an
 *   optional operation this is reported as a compile error (see the message
 *   "different signature than the base class declaration").
 *
 * - Optional operations must be overridden as public members, or the derived
 *   class must declare 'friend sundials::NVectorBase<MyVec>;'. Marking the
 *   derived class final lets the compiler call overrides directly when the
 *   static type is known.
 *
 * - Each operation name must be a single function in the derived class, the
 *   override. Do not declare overloaded or template helpers with an operation's
 *   name (e.g., Print(std::ostream&)); this is a compile error at the
 *   operation's line in sunInstallOps.
 *
 * - Clone must allocate the new object with new (N_VDestroy deletes it), and
 *   the destructor must be public.
 *
 * - std::vector<MyVec> growth keeps handles obtained earlier valid only if
 *   MyVec's move constructor is noexcept (so all data members must be
 *   nothrow-movable). Otherwise std::vector copies on reallocation, the copies
 *   get new handles, and handles obtained earlier dangle. Check it in your own
 *   code with static_assert(std::is_nothrow_move_constructible_v<MyVec>);
 *
 * - Rule of five. A derived class that declares a destructor or any copy member
 *   loses its implicit moves, and moves then silently copy. It must explicitly
 *   default its move constructor and move assignment. A hand-written copy
 *   constructor must call NVectorBase(other), otherwise the copy has no N_Vector
 *   (GCC -Wextra and clang-tidy bugprone-copy-constructor-init warn). Prefer
 *   defaulted copy and move members.
 *
 * -----------------------------------------------------------------------------
 * OPERATIONS
 * -----------------------------------------------------------------------------
 *
 * Every operation falls in one of three categories. The declarations in the
 * class below are the reference for the exact signatures. The view types
 * (ConstVectors, Reals, ...) are described in ARRAY ARGUMENTS below.
 *
 * 1. Required. Public and pure virtual: a derived class that does not override
 *    them cannot be instantiated. Always installed.
 *      Clone, GetLength, LinearSum, Const, Scale, DotProd, MaxNorm
 *
 * 2. Real default. Public and virtual, always installed (N_VGetVectorID calls
 *    it without a NULL check).
 *      GetVectorID returns SUNDIALS_NVEC_CUSTOM. The default must not throw.
 *
 * 3. Optional. Protected and virtual, installed only if the derived class
 *    overrides the operation. The default throws std::logic_error. An operation
 *    that is not overridden stays NULL in the ops table, so the SUNDIALS C
 *    fallbacks (e.g., the fused and vector array operations) or the integrator
 *    checks at initialization ("A required vector operation is not
 *    implemented.") apply.
 *      GetCommunicator, GetLocalLength, GetArrayPointer, SetArrayPointer,
 *      GetDeviceArrayPointer, SetDeviceArrayPointer
 *      Prod, Div, Abs, Inv, AddConst, Compare, WrmsNorm, WrmsNormMask, Min,
 *      WL2Norm, L1Norm, MinQuotient, InvTest, ConstrMask
 *      LinearCombination, ScaleAddMulti, DotProdMulti
 *      LinearSumVectorArray, ScaleVectorArray, ConstVectorArray,
 *      WrmsNormVectorArray, WrmsNormMaskVectorArray, ScaleAddMultiVectorArray,
 *      LinearCombinationVectorArray
 *      DotProdLocal, MaxNormLocal, MinLocal, L1NormLocal, InvTestLocal,
 *      ConstrMaskLocal, MinQuotientLocal, WSqrSumLocal, WSqrSumMaskLocal
 *      DotProdMultiLocal, DotProdMultiAllReduce
 *      BufSize, BufPack, BufUnpack
 *      Print, PrintFile
 *
 * The destroy operation is always installed too, but it is not an operation a
 * derived class can override.
 *
 * Effects of protected optional operations with a throwing default:
 *
 * - Through SUNDIALS (the C API): a non-overridden operation has no trampoline,
 *   so C code never reaches the default.
 * - Through user C++ code: myvec.Div(x, y) compiles only if MyVec overrides Div
 *   as public.
 * - A derived class member that calls an operation it does not override reaches
 *   the default and throws std::logic_error. Inside a trampoline, that is
 *   reported as SUN_ERR_OP_FAIL (see ERRORS AND EXCEPTIONS).
 * - The default of an optional operation can be called as notImplemented("Op")
 *   from a derived class that wants to install an operation but fail at run
 *   time.
 *
 * GetVectorID() may be overridden. Do not return the ID of a built-in vector
 * type (e.g., SUNDIALS_NVEC_SERIAL): SUNDIALS modules may use the ID to decide
 * how to interpret the content pointer.
 *
 * Optional operations needed by the SUNDIALS packages, beyond the 7 required
 * operations (a missing one is reported at initialization or when the feature
 * is used):
 *
 *   CVODE(S)   Prod, Div, Abs, Inv, AddConst, WrmsNorm; plus Min with vector
 *              tolerances or atolmin0
 *   ARKODE     Div, Abs, Inv, WrmsNorm; plus Min and AddConst depending on the
 *              tolerance type
 *   IDA(S)     Prod, Abs, Inv, AddConst, WrmsNorm, Min; plus WrmsNormMask with
 *              suppressalg
 *   KINSOL     Prod, Div, Abs, Inv, Min, WL2Norm; plus L1Norm for the
 *              difference quotient Jacobian-times-vector
 *   Constraints (any package)
 *              Compare, ConstrMask, MinQuotient, Div, MaxNorm; plus Prod in
 *              ARKODE
 *   CVDiag     Compare, InvTest
 *   Used as a ManyVector, MPIManyVector or MPIPlusX subvector: every operation
 *              the package lists, plus GetLocalLength
 *
 * -----------------------------------------------------------------------------
 * SUBCLASSES OF THE DERIVED CLASS
 * -----------------------------------------------------------------------------
 *
 * final is not required. The rule is that an optional operation is installed
 * only if the class Derived in NVectorBase<Derived> itself declares it. A subclass of Derived
 * can override the operations that Derived installs, and they dispatch virtually.
 * To let a subclass add an optional operation, Derived must override that operation
 * too (with a real implementation, or one that calls notImplemented):
 *
 *   class VecA : public sundials::NVectorBase<VecA> { ... Div ... };
 *   class VecB : public VecA
 *   {
 *   public:
 *     void Prod(const VecA& x, const VecA& y) override; // note: VecA
 *   };
 *
 * Here a VecB object is built through NVectorBase<VecA>, which inspects VecA.
 * VecA doesn't declare Prod, so nvprod stays NULL for VecB objects as well and
 * SUNDIALS never calls VecB::Prod; the integrator rejects the vector at
 * initialization. Declaring Prod in VecA fixes this. Marking a class final
 * turns this silent case into a compile error.
 *
 * -----------------------------------------------------------------------------
 * BASE CLASS DESIGN AND TRAMPOLINES
 * -----------------------------------------------------------------------------
 *
 * 1. CRTP (Curiously Recurring Template Pattern)
 *
 *    The derived class passes itself as a template parameter to the base class:
 *
 *    class MyVec : public NVectorBase<MyVec> { ... };
 *
 *    This lets the base class declare operations in terms of the derived type,
 *    e.g., virtual void LinearSum(sunrealtype, const MyVec&, sunrealtype,
 *    const MyVec&), and cast the N_Vector void* content pointer to MyVec*. The
 *    content pointer stores the address of the impl::ObjectBase subobject (the
 *    handle is created in the base class constructor, before the MyVec object
 *    exists) and FromNVector() performs the downcast to MyVec* when an
 *    operation is called.
 *
 * 2. VIRTUAL OPERATIONS
 *
 *    Each operation is a virtual function, so a call through the base class
 *    reaches the derived class override even when the derived class is
 *    further subclassed.
 *
 * 3. STATIC TRAMPOLINES (function pointer bridge)
 *
 *    The N_Vector ops table is a struct of C function pointers, e.g.:
 *
 *    void (*nvscale)(sunrealtype, N_Vector, N_Vector);
 *
 *    A C function pointer cannot point to a C++ non-static member function, so
 *    we define static member functions ("trampolines") with the C signature,
 *    extract the C++ object, and call the corresponding virtual function. For
 *    example, the ops table stores
 *
 *    ops->nvscale = &nvScale;
 *
 *    where the static function is (with the exception guard described in
 *    ERRORS AND EXCEPTIONS below)
 *
 *    static void nvScale(sunrealtype c, N_Vector x, N_Vector z) noexcept
 *    {
 *      impl::guardVoid(z->sunctx, "N_VScale", __FILE__, __LINE__,
 *                      [&] { base(z)->Scale(c, *FromNVectorConst(x)); });
 *    }
 *
 *    base(z) is FromNVector(z) converted to NVectorBase*, so the call is
 *    dispatched through the base class. The convention is to use the output
 *    vector as "this" and the input vectors as const MyVec& arguments.
 *    Operations without an output vector (norms, dot products, etc.) are const
 *    methods called on the first input vector.
 *
 * -----------------------------------------------------------------------------
 * OWNERSHIP AND LIFETIME
 * -----------------------------------------------------------------------------
 *
 * Every NVectorBase object constructed with a SUNContext owns one N_Vector
 * handle whose content pointer points back at the object. The handle is
 * created by the constructor and freed by the destructor. There are two kinds
 * of objects:
 *
 * - User-owned objects are created by user code (on the stack, with new,
 *   in containers, ...) and destroyed by user code. The handle can be passed
 *   to any SUNDIALS function. N_VDestroy on such a handle does nothing (it is
 *   silent), and the object's destructor frees the handle.
 *
 * - Clones are created by N_VClone, which calls Clone() and marks the new
 *   object as owned by its handle. N_VDestroy on the handle deletes the object
 *   (which frees the handle), so Clone() must allocate with new. Only clones
 *   are deleted by N_VDestroy.
 *
 * Copy, move and assignment follow these rules (the semantics live in
 * impl::ObjectBase, see sundials_object_base.hpp):
 *
 * - Copying an object without a handle gives an object without a handle.
 *   Otherwise the copy gets a new handle from the source's SUNContext (this may
 *   throw), and the source is unchanged.
 * - Moving from a user-owned object with a handle takes the handle and
 *   repoints its content to the new object, so an N_Vector obtained from the
 *   source now refers to the target. The source is left without a handle.
 * - Moving from a clone allocates a new handle for the new object, because the
 *   clone's handle belongs to SUNDIALS and the clone keeps it. If that
 *   allocation fails, the error is reported through the source's SUNContext
 *   (SUN_ERR_MEM_FAIL) and the new object has no handle; its data is still
 *   moved.
 * - Assignment (copy or move) never removes the target's handle. A target
 *   with a handle keeps it and only the data is assigned, even from a source
 *   without a handle (the target then is a live vector holding that data, e.g.,
 *   of length 0). A target without a handle gets one as in the corresponding
 *   constructor.
 * - The clone flag is never copied, moved or transferred.
 * - The default constructor and both moves are noexcept; no path terminates.
 *
 * Consequences:
 *
 * - If MyVec's move constructor is noexcept (all data members nothrow-movable;
 *   check with static_assert(std::is_nothrow_move_constructible_v<MyVec>)),
 *   std::vector<MyVec> growth moves the elements and keeps handles obtained
 *   earlier valid. Otherwise std::vector copies on reallocation and those
 *   handles dangle. std::swap exchanges handles along with
 *   the data (std::swap with a clone exchanges data only).
 * - *MyVec::FromNVector(ydot) = MakeRhs(); and
 *   *MyVec::FromNVector(ydot) = MyVec{}; are both safe. After the second, ydot
 *   is live with the empty object's data, so later operations may see
 *   mismatched lengths; that is the user's responsibility.
 * - Moving out of an N_VClone result is safe.
 * - Assignment keeps the target's SUNContext, even if the source was created
 *   with another one.
 *
 * Objects without a handle: a default-constructed object, a moved-from
 * non-clone, and (after an out-of-memory error) a moved-to object have no
 * handle. HasNVector() returns false and GetNVector() returns NULL (a debug
 * check failure). Passing NULL to any N_V function crashes, so test
 * HasNVector() first. The only valid operations on a moved-from object are
 * destruction and assignment.
 *
 * Reusing allocations: there is no hook for recycling cloned objects. Expensive
 * resources (e.g., device memory or a distributed vector) can be pooled by
 * holding them in a data member that returns them to a pool when destroyed
 * (e.g., a smart pointer with a custom deleter); N_VDestroy -> delete runs the
 * derived class destructor, which releases the member. The memory of the
 * derived object itself can be recycled with class-specific operator new and
 * operator delete, which Clone() (new) and N_VDestroy (delete) use. In both
 * cases the N_Vector handle is still created and freed for every clone.
 *
 * -----------------------------------------------------------------------------
 * ARRAY ARGUMENTS
 * -----------------------------------------------------------------------------
 *
 * Array arguments are passed as non-owning views that wrap the arrays SUNDIALS
 * passes in without allocating or copying. NVectorBase<Derived> provides the aliases
 *
 *   ConstVectors   = VectorSpan<const Derived>    N_Vector*  inputs,  X[j] -> const Derived&
 *   Vectors        = VectorSpan<Derived>          N_Vector*  outputs, Z[j] -> Derived&
 *   ConstVectors2D = VectorSpan2D<const Derived>  N_Vector** inputs,  Y[j][i]
 *   Vectors2D      = VectorSpan2D<Derived>        N_Vector** outputs, Z[j][i]
 *   ConstReals     = ScalarSpan<const sunrealtype>  scalar coefficient inputs
 *   Reals          = ScalarSpan<sunrealtype>        scalar result outputs
 *
 * All views provide size() and operator[]; the 1D views also support
 * range-based for loops. For 2D views, the outer index is the row (e.g., nsum)
 * and the inner index is the column (e.g., nvec), matching the SUNDIALS
 * N_Vector** layout, and rows() and cols() give the dimensions. data() returns
 * the underlying pointer, e.g., to detect in-place calls (Y.data() == Z.data()).
 *
 * The views are valid only for the duration of the call. ScalarSpan could be
 * replaced by std::span when C++20 can be required.
 *
 * -----------------------------------------------------------------------------
 * ALIASING
 * -----------------------------------------------------------------------------
 *
 * SUNDIALS allows output vectors to be the same objects as input vectors, so a
 * const Derived& argument (or an element of a ConstVectors view) may refer to *this
 * or to an element of an output view. For example:
 *
 * - LinearSum: z may be x and/or y. The SUNDIALS C fallback for
 *   LinearCombination (sundials_nvector.c) calls nvlinearsum(c[i], X[i], 1, z,
 *   z), so z is also an input there.
 * - LinearCombination: X[0] may be z (commonly with c[0] == 1).
 * - ScaleAddMulti and ScaleAddMultiVectorArray: Z may be Y (in-place).
 * - ScaleVectorArray: Z may be X (in-place).
 * - LinearCombinationVectorArray: X[0] may be Z.
 *
 * Element-wise loops that read all inputs at index i before writing the output
 * at index i handle aliasing automatically. Implementations that make multiple
 * passes, resize, or use deferred evaluation must account for it.
 *
 * -----------------------------------------------------------------------------
 * ERRORS AND EXCEPTIONS
 * -----------------------------------------------------------------------------
 *
 * Exceptions must not propagate into SUNDIALS C code. Every trampoline is
 * noexcept and catches any exception thrown by an operation, including the
 * std::logic_error thrown by an optional operation that is not overridden.
 * The exception is reported through the SUNContext of the dispatch vector
 * (this records the last error and runs the context's error handlers) with the
 * name of the C API as the function (e.g., "N_VLinearSum"), the what() message,
 * and the code
 *
 * - SUN_ERR_MEM_FAIL for std::bad_alloc;
 * - SUN_ERR_OP_FAIL for any other exception ("unknown C++ exception" if it is
 *   not a std::exception).
 *
 * Clone returning NULL, or an object without an N_Vector, is also reported as
 * SUN_ERR_OP_FAIL (N_VCloneVectorArray detects failure only through the last
 * error). After reporting an error, the trampoline returns:
 *
 *   SUNErrCode ops (fused, vector array, single buffer, XBraid)  the error code
 *   sunbooleantype ops (InvTest, ConstrMask, InvTestLocal,
 *     ConstrMaskLocal)                                           SUNFALSE
 *   void ops (element-wise, print)                               nothing
 *   GetVectorID                                                  SUNDIALS_NVEC_CUSTOM
 *   GetCommunicator                                              SUN_COMM_NULL
 *   norms, dot products, Min, lengths                            0
 *   array pointers, Clone                                        NULL
 *
 * A placeholder is returned only because the C signature requires a value. It
 * is not an error indicator, and placeholders come only from exceptions: an
 * integrator continues with the placeholder after the error is logged at the
 * default logging level. For a hard stop, push SUNAbortErrHandlerFn onto the
 * context's error handler stack. A handler that throws inside a noexcept
 * trampoline (or a derived class destructor that throws inside N_VDestroy)
 * ends in std::terminate.
 *
 * The what() message is reported verbatim: a % in it is not a format directive
 * and long messages are not truncated.
 *
 * A derived class reports an error by throwing an exception (the trampoline
 * reports it), or, for operations returning SUNErrCode, by returning an error
 * code. SUNDIALS checks the last error with the SUNCheckLastErr* macros when
 * SUNDIALS_ENABLE_ERROR_CHECKS is defined.
 *
 * To report an error without throwing, call
 *
 *   SUNContext_HandleError(GetContext(), code, __LINE__, __func__, __FILE__,
 *                          "msg");
 *
 * and then return normally (void and value operations) or return the code
 * (SUNErrCode operations). Throwing stays the primary path. Report or throw,
 * not both, or the error is reported twice. A handler that throws inside a
 * noexcept trampoline ends in std::terminate.
 *
 * -----------------------------------------------------------------------------
 * DEBUG CHECKS
 * -----------------------------------------------------------------------------
 *
 * When SUNDIALS_DEBUG is defined, the following misuses print a message to
 * stderr and call std::abort(). Otherwise they are not checked. SUNDIALS_DEBUG
 * must be defined when compiling the code that includes this header (e.g.,
 * -DSUNDIALS_DEBUG); the SUNDIALS_DEBUG CMake option only applies to the
 * SUNDIALS library sources. It must be defined consistently across all
 * translation units of a program: mixing translation units compiled with and
 * without it gives two definitions of the same inline functions (an ODR
 * violation).
 *
 * - GetNVector() on an object without a handle.
 * - FromNVector() on a NULL or empty N_Vector, or on an N_Vector that was not
 *   created by NVectorBase<Derived> (e.g., mixing two custom vector types). The
 *   check compares the handle's nvdestroy pointer with this class's destroy
 *   trampoline. Linkers that fold identical functions (e.g., MSVC /OPT:ICF or
 *   lld --icf=all) may merge the trampolines of different Derived types, which
 *   can hide a mismatch but never reports a false one.
 * - Out of range indices in the array views.
 *
 * -----------------------------------------------------------------------------
 * LIMITATIONS
 * -----------------------------------------------------------------------------
 *
 * - N_VCloneEmpty and the nvspace operation are not provided (they stay NULL).
 *   Because N_VCloneEmpty is not provided, any path that clones a ManyVector or
 *   MPIManyVector "empty" with an NVectorBase subvector crashes. The in-tree
 *   callers are the difference quotient Jacobian code, which also needs array
 *   pointers that ManyVector does not provide, so this path is unlikely.
 * - Unchecked uses of optional operations. ManyVector, MPIManyVector and
 *   MPIPlusX populate a full ops table regardless of what their subvectors
 *   provide and forward to the subvector operations without checking, so a
 *   subvector that lacks an operation the package needs is a NULL call (see the
 *   table in OPERATIONS; GetLocalLength and N_VCloneEmpty are also called on
 *   subvectors unchecked). The same holds for N_VL1Norm in the dominant
 *   eigenvalue estimators, and for the Fortran bindings, which call the N_V
 *   functions directly. The iterative linear solvers (SPGMR, SPFGMR, SPBCGS,
 *   SPTFQMR) use Prod and Div only when scaling vectors are given, which
 *   in-tree happens only from KINSOL, whose initialization checks for both.
 *   However, in builds with error checks (SUNDIALS_ENABLE_ERROR_CHECKS), these
 *   solvers require Prod and Div at construction (e.g., sunlinsol_spgmr.c:87),
 *   even when no scaling vectors are used.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_SUNDIALS_NVECTOR_BASE_HPP
#define SUNDIALS_SUNDIALS_NVECTOR_BASE_HPP

#if defined(_MSVC_LANG)
#if _MSVC_LANG < 201703L
#error "sundials_nvector_base.hpp requires C++17 or later"
#endif
#elif !defined(__cplusplus) || __cplusplus < 201703L
#error "sundials_nvector_base.hpp requires C++17 or later"
#endif

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <sundials/sundials_object_base.hpp>

#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

// Debug check used throughout this header (see DEBUG CHECKS above). Undefined
// at the end of this header.
#if defined(SUNDIALS_DEBUG)
#define SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(cond, msg)                        \
  do {                                                                      \
    if (!(cond))                                                            \
    {                                                                       \
      std::fprintf(stderr, "%s:%d: %s: %s\n", __FILE__, __LINE__, __func__, \
                   msg);                                                    \
      std::abort();                                                         \
    }                                                                       \
  }                                                                         \
  while (0)
#else
#define SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(cond, msg) ((void)0)
#endif

namespace sundials {

template<class Derived>
class NVectorBase;

// =============================================================================
// Non-owning views of the arrays SUNDIALS passes to fused and vector array
// operations. See ARRAY ARGUMENTS above.
//
// Style note for code review: the view member functions use standard library
// (lowercase snake_case) names rather than the SUNDIALS PascalCase convention
// for public C++ functions. begin() and end() are required by range-based for
// loops, and size(), data(), empty(), and operator[] match std::span and work
// with generic code such as std::size and std::data. rows() and cols() follow
// the same style for consistency within the view types.
// =============================================================================

// View of a scalar array. ScalarType is sunrealtype or const sunrealtype.
// With C++20 this could be replaced by std::span<ScalarType>.
template<class ScalarType>
class ScalarSpan
{
public:
  using element_type = ScalarType;
  using value_type   = std::remove_cv_t<ScalarType>;
  using size_type    = std::size_t;
  using iterator     = ScalarType*;

  constexpr ScalarSpan() noexcept = default;

  constexpr ScalarSpan(ScalarType* data, size_type size) noexcept
    : data_(data), size_(size)
  {}

  constexpr ScalarType* data() const noexcept { return data_; }

  constexpr size_type size() const noexcept { return size_; }

  constexpr bool empty() const noexcept { return size_ == 0; }

  constexpr ScalarType& operator[](size_type i) const noexcept
  {
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(i < size_,
                                      "ScalarSpan index out of range");
    return data_[i];
  }

  constexpr iterator begin() const noexcept { return data_; }

  constexpr iterator end() const noexcept { return data_ + size_; }

private:
  ScalarType* data_        = nullptr;
  size_type size_ = 0;
};

// View of an N_Vector* array presented as derived objects. VectorType is Derived (for
// outputs) or const Derived (for inputs).
template<class VectorType>
class VectorSpan
{
public:
  using element_type = VectorType;
  using value_type   = std::remove_cv_t<VectorType>;
  using size_type    = std::size_t;

  class iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = std::remove_cv_t<VectorType>;
    using difference_type   = std::ptrdiff_t;
    using pointer           = VectorType*;
    using reference         = VectorType&;

    iterator() noexcept = default;

    explicit iterator(N_Vector* ptr) noexcept : ptr_(ptr) {}

    reference operator*() const noexcept { return VectorSpan::element(*ptr_); }

    pointer operator->() const noexcept { return &VectorSpan::element(*ptr_); }

    iterator& operator++() noexcept
    {
      ++ptr_;
      return *this;
    }

    iterator operator++(int) noexcept
    {
      iterator tmp = *this;
      ++ptr_;
      return tmp;
    }

    friend bool operator==(const iterator& lhs, const iterator& rhs) noexcept
    {
      return lhs.ptr_ == rhs.ptr_;
    }

    friend bool operator!=(const iterator& lhs, const iterator& rhs) noexcept
    {
      return lhs.ptr_ != rhs.ptr_;
    }

  private:
    N_Vector* ptr_ = nullptr;
  };

  VectorSpan() noexcept = default;

  VectorSpan(N_Vector* data, size_type size) noexcept : data_(data), size_(size)
  {}

  // Underlying N_Vector array, e.g., to pass on to other SUNDIALS functions or
  // to detect in-place calls.
  N_Vector* data() const noexcept { return data_; }

  size_type size() const noexcept { return size_; }

  bool empty() const noexcept { return size_ == 0; }

  VectorType& operator[](size_type i) const noexcept
  {
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(i < size_,
                                      "VectorSpan index out of range");
    return element(data_[i]);
  }

  iterator begin() const noexcept { return iterator(data_); }

  iterator end() const noexcept { return iterator(data_ + size_); }

private:
  static VectorType& element(N_Vector vec) noexcept
  {
    return *NVectorBase<std::remove_cv_t<VectorType>>::FromNVector(vec);
  }

  N_Vector* data_ = nullptr;
  size_type size_ = 0;
};

// View of an N_Vector** array (rows x cols) presented as derived objects. VectorType is
// Derived (for outputs) or const Derived (for inputs).
template<class VectorType>
class VectorSpan2D
{
public:
  using size_type = std::size_t;

  VectorSpan2D() noexcept = default;

  VectorSpan2D(N_Vector** data, size_type rows, size_type cols) noexcept
    : data_(data), rows_(rows), cols_(cols)
  {}

  N_Vector** data() const noexcept { return data_; }

  size_type rows() const noexcept { return rows_; }

  size_type cols() const noexcept { return cols_; }

  // The size is the number of rows, i.e., the number of 1D views in the outer
  // index (use cols() for the length of each row)
  size_type size() const noexcept { return rows_; }

  bool empty() const noexcept { return rows_ == 0; }

  // Row j as a 1D view with cols() elements
  VectorSpan<VectorType> operator[](size_type j) const noexcept
  {
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(j < rows_,
                                      "VectorSpan2D row index out of range");
    return VectorSpan<VectorType>(data_[j], cols_);
  }

private:
  N_Vector** data_ = nullptr;
  size_type rows_  = 0;
  size_type cols_  = 0;
};

namespace impl {

// Per-type details for ObjectBase (see sundials_object_base.hpp)
struct NVectorTraits
{
  using Struct = std::remove_pointer_t<N_Vector>;

  static N_Vector NewEmpty(SUNContext ctx) { return N_VNewEmpty(ctx); }

  static void FreeEmpty(N_Vector vec) { N_VFreeEmpty(vec); }

  static constexpr const char* CLASS_NAME = "sundials::NVectorBase"; // for error messages
};

} // namespace impl

// =============================================================================
// NVectorBase<Derived>
//
// CRTP base class with virtual operations. Owns the N_Vector handle for the
// derived object.
// =============================================================================

template<class Derived>
class NVectorBase
  // The third argument is the class itself. ObjectBase calls the installer as
  // NVectorBase<Derived>::sunInstallOps (a qualified call), which is what stops
  // a member of Derived with the same name from replacing it.
  : public impl::ObjectBase<Derived, impl::NVectorTraits, NVectorBase<Derived>>
{
private:
  using Base = impl::ObjectBase<Derived, impl::NVectorTraits, NVectorBase>;
  // ObjectBase befriends this class (friend TypeBase) so that this class can
  // use its private members, and this class befriends ObjectBase so that
  // ObjectBase can call the private sunInstallOps. Keep both.
  friend Base;

public:
  // ---------------------------------------------------------------------------
  // Array argument views
  // ---------------------------------------------------------------------------

  using ConstVectors   = VectorSpan<const Derived>;
  using Vectors        = VectorSpan<Derived>;
  using ConstVectors2D = VectorSpan2D<const Derived>;
  using Vectors2D      = VectorSpan2D<Derived>;
  using ConstReals     = ScalarSpan<const sunrealtype>;
  using Reals          = ScalarSpan<sunrealtype>;

  // ---------------------------------------------------------------------------
  // Construction, copy, move and destruction
  //
  // The default constructor creates an object without an N_Vector handle, e.g.,
  // for containers or two-phase initialization. The SUNContext constructor
  // creates the handle and throws std::invalid_argument if ctx is NULL and
  // std::bad_alloc if the handle cannot be allocated.
  //
  // The copy and move members are protected (see below), so a derived object
  // cannot be sliced through a reference to the base class.
  // ---------------------------------------------------------------------------

  NVectorBase() noexcept = default;

  explicit NVectorBase(SUNContext ctx) : Base(ctx) {}

  virtual ~NVectorBase() = default;

protected:
  // Copy and move are protected to prevent slicing (C.67 in the C++ Core
  // Guidelines: a polymorphic class should suppress public copy and move).
  // They are defaulted, because the virtual destructor would otherwise suppress
  // the implicit moves, and the handle semantics live in impl::ObjectBase.
  // Derived classes use them for their own defaulted copy and move members.
  NVectorBase(const NVectorBase&) = default;

  NVectorBase(NVectorBase&&) noexcept = default;

  NVectorBase& operator=(const NVectorBase&) = default;

  NVectorBase& operator=(NVectorBase&&) noexcept = default;

public:

  // ---------------------------------------------------------------------------
  // Handle access
  //
  // The N_Vector handle can be passed directly to any C function expecting an
  // N_Vector. These never allocate.
  //
  // Both are deleted for rvalues (temporaries): the handle is freed when the
  // temporary is destroyed at the end of the full expression, so passing it to
  // a SUNDIALS function that keeps the handle (e.g., N_VNew_ManyVector keeps
  // its subvectors) would leave SUNDIALS with a dangling pointer.
  //
  // Note that the C API has no notion of const, so the handle of a const
  // object allows SUNDIALS to modify the object.
  //
  // Returns NULL for an object without a handle (a debug check failure).
  // ---------------------------------------------------------------------------

  N_Vector GetNVector() const& noexcept
  {
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(this->objHasHandle(),
                                      "NVectorBase: vector has no N_Vector "
                                      "(default-constructed or moved-from)");
    return this->objHandle();
  }

  N_Vector GetNVector() const&& = delete;

  operator N_Vector() const& noexcept { return GetNVector(); }

  operator N_Vector() const&& = delete;

  // True if the object has an N_Vector handle, i.e., it is not
  // default-constructed or moved-from. Use this rather than testing the
  // N_Vector conversion (e.g., if (vec)), which calls GetNVector() and fails the
  // debug check for an object without a handle.
  bool HasNVector() const noexcept { return this->objHasHandle(); }

  // The SUNContext of the handle (NULL if there is no handle), e.g., to create
  // other SUNDIALS objects with the same context.
  SUNContext GetContext() const noexcept { return this->objContext(); }

  // ---------------------------------------------------------------------------
  // Derived object access
  //
  // Recover a typed Derived pointer from the void* content field of an
  // N_Vector. These are public so that user code outside the class hierarchy
  // can use them e.g., in callback functions that receive N_Vector arguments
  // and need to access the underlying C++ data (ODE right-hand sides, Jacobian
  // evaluations, preconditioner solves, etc.):
  //
  // int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
  // {
  //   const auto* yv = MyVec::FromNVectorConst(y);
  //   auto* ydotv    = MyVec::FromNVector(ydot);
  //   // ... work directly with MyVec members ...
  // }
  //
  // When SUNDIALS_DEBUG is defined, these check that the N_Vector was created
  // by NVectorBase<Derived> (see DEBUG CHECKS above).
  // ---------------------------------------------------------------------------

  // Get a mutable Derived pointer from an N_Vector.
  static Derived* FromNVector(N_Vector vec) noexcept
  {
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(vec != nullptr && vec->content != nullptr &&
                                        vec->ops != nullptr,
                                      "NVectorBase: NULL or empty N_Vector");
    SUNDIALS_NVECTOR_BASE_DEBUG_CHECK(vec->ops->nvdestroy == &nvDestroy,
                                      "NVectorBase: N_Vector was not created "
                                      "by this vector type");
    return Base::objFromHandle(vec);
  }

  // Get a const Derived pointer from an N_Vector.
  static const Derived* FromNVectorConst(N_Vector vec) noexcept
  {
    return FromNVector(vec);
  }

  // ===========================================================================
  // Operations
  //
  // The comments give the operation computed; "this" is the output vector z
  // (or the first input for operations without an output vector).
  // ===========================================================================

  // ---------------------------------------------------------------------------
  // Category 1: required operations
  // ---------------------------------------------------------------------------

  // Used by N_VClone: allocate (with new) a vector with the same structure.
  // The values of the new vector do not need to be initialized.
  virtual Derived* Clone(SUNContext ctx) const = 0;

  // Global length
  virtual sunindextype GetLength() const = 0;

  // z = a*x + b*y
  virtual void LinearSum(sunrealtype a, const Derived& x, sunrealtype b,
                         const Derived& y) = 0;

  // z_i = c
  virtual void Const(sunrealtype c) = 0;

  // z = c*x
  virtual void Scale(sunrealtype c, const Derived& x) = 0;

  // Returns x . y
  virtual sunrealtype DotProd(const Derived& y) const = 0;

  // Returns max_i |x_i|
  virtual sunrealtype MaxNorm() const = 0;

  // ---------------------------------------------------------------------------
  // Category 2: operation with a real default
  // ---------------------------------------------------------------------------

  virtual N_Vector_ID GetVectorID() const { return SUNDIALS_NVEC_CUSTOM; }

protected:
  // ===========================================================================
  // Category 3: optional operations
  //
  // Installed only if Derived overrides them. The defaults throw
  // std::logic_error through notImplemented.
  // ===========================================================================

  // ---------------------------------------------------------------------------
  // Utility operations
  // ---------------------------------------------------------------------------

  virtual SUNComm GetCommunicator() const
  {
    notImplemented("GetCommunicator");
  }

  virtual sunindextype GetLocalLength() const
  {
    notImplemented("GetLocalLength");
  }

  virtual sunrealtype* GetArrayPointer()
  {
    notImplemented("GetArrayPointer");
  }

  virtual void SetArrayPointer(sunrealtype* /* data */)
  {
    notImplemented("SetArrayPointer");
  }

  virtual sunrealtype* GetDeviceArrayPointer()
  {
    notImplemented("GetDeviceArrayPointer");
  }

  virtual void SetDeviceArrayPointer(sunrealtype* /* data */)
  {
    notImplemented("SetDeviceArrayPointer");
  }

  // ---------------------------------------------------------------------------
  // Standard operations
  // ---------------------------------------------------------------------------

  virtual void Prod(const Derived& /* x */, const Derived& /* y */)
  {
    notImplemented("Prod");
  }

  virtual void Div(const Derived& /* x */, const Derived& /* y */)
  {
    notImplemented("Div");
  }

  virtual void Abs(const Derived& /* x */)
  {
    notImplemented("Abs");
  }

  virtual void Inv(const Derived& /* x */)
  {
    notImplemented("Inv");
  }

  virtual void AddConst(const Derived& /* x */, sunrealtype /* b */)
  {
    notImplemented("AddConst");
  }

  virtual void Compare(sunrealtype /* c */, const Derived& /* x */)
  {
    notImplemented("Compare");
  }

  virtual sunrealtype WrmsNorm(const Derived& /* w */) const
  {
    notImplemented("WrmsNorm");
  }

  virtual sunrealtype WrmsNormMask(const Derived& /* w */, const Derived& /* id */) const
  {
    notImplemented("WrmsNormMask");
  }

  virtual sunrealtype Min() const
  {
    notImplemented("Min");
  }

  virtual sunrealtype WL2Norm(const Derived& /* w */) const
  {
    notImplemented("WL2Norm");
  }

  virtual sunrealtype L1Norm() const
  {
    notImplemented("L1Norm");
  }

  virtual sunrealtype MinQuotient(const Derived& /* denom */) const
  {
    notImplemented("MinQuotient");
  }

  virtual sunbooleantype InvTest(const Derived& /* x */)
  {
    notImplemented("InvTest");
  }

  virtual sunbooleantype ConstrMask(const Derived& /* c */, const Derived& /* x */)
  {
    notImplemented("ConstrMask");
  }

  // ---------------------------------------------------------------------------
  // Fused operations
  // ---------------------------------------------------------------------------

  virtual SUNErrCode LinearCombination(ConstReals /* c */, ConstVectors /* X */)
  {
    notImplemented("LinearCombination");
  }

  virtual SUNErrCode ScaleAddMulti(ConstReals /* a */, ConstVectors /* Y */, Vectors /* Z */) const
  {
    notImplemented("ScaleAddMulti");
  }

  virtual SUNErrCode DotProdMulti(ConstVectors /* Y */, Reals /* dotprods */) const
  {
    notImplemented("DotProdMulti");
  }

  // ---------------------------------------------------------------------------
  // Vector array operations
  // ---------------------------------------------------------------------------

  virtual SUNErrCode LinearSumVectorArray(sunrealtype /* a */, ConstVectors /* X */, sunrealtype /* b */, ConstVectors /* Y */, Vectors /* Z */)
  {
    notImplemented("LinearSumVectorArray");
  }

  virtual SUNErrCode ScaleVectorArray(ConstReals /* c */, ConstVectors /* X */, Vectors /* Z */)
  {
    notImplemented("ScaleVectorArray");
  }

  virtual SUNErrCode ConstVectorArray(sunrealtype /* c */, Vectors /* Z */)
  {
    notImplemented("ConstVectorArray");
  }

  virtual SUNErrCode WrmsNormVectorArray(ConstVectors /* X */, ConstVectors /* W */, Reals /* nrm */) const
  {
    notImplemented("WrmsNormVectorArray");
  }

  virtual SUNErrCode WrmsNormMaskVectorArray(ConstVectors /* X */, ConstVectors /* W */, const Derived& /* id */, Reals /* nrm */) const
  {
    notImplemented("WrmsNormMaskVectorArray");
  }

  virtual SUNErrCode ScaleAddMultiVectorArray(ConstReals /* a */, ConstVectors /* X */, ConstVectors2D /* Y */, Vectors2D /* Z */) const
  {
    notImplemented("ScaleAddMultiVectorArray");
  }

  virtual SUNErrCode LinearCombinationVectorArray(ConstReals /* c */, ConstVectors2D /* X */, Vectors /* Z */)
  {
    notImplemented("LinearCombinationVectorArray");
  }

  // ---------------------------------------------------------------------------
  // Local reduction operations
  // ---------------------------------------------------------------------------

  virtual sunrealtype DotProdLocal(const Derived& /* y */) const
  {
    notImplemented("DotProdLocal");
  }

  virtual sunrealtype MaxNormLocal() const
  {
    notImplemented("MaxNormLocal");
  }

  virtual sunrealtype MinLocal() const
  {
    notImplemented("MinLocal");
  }

  virtual sunrealtype L1NormLocal() const
  {
    notImplemented("L1NormLocal");
  }

  virtual sunbooleantype InvTestLocal(const Derived& /* x */)
  {
    notImplemented("InvTestLocal");
  }

  virtual sunbooleantype ConstrMaskLocal(const Derived& /* c */, const Derived& /* x */)
  {
    notImplemented("ConstrMaskLocal");
  }

  virtual sunrealtype MinQuotientLocal(const Derived& /* denom */) const
  {
    notImplemented("MinQuotientLocal");
  }

  virtual sunrealtype WSqrSumLocal(const Derived& /* w */) const
  {
    notImplemented("WSqrSumLocal");
  }

  virtual sunrealtype WSqrSumMaskLocal(const Derived& /* w */, const Derived& /* id */) const
  {
    notImplemented("WSqrSumMaskLocal");
  }

  // ---------------------------------------------------------------------------
  // Single buffer reduction operations
  // ---------------------------------------------------------------------------

  virtual SUNErrCode DotProdMultiLocal(ConstVectors /* Y */, Reals /* dotprods */) const
  {
    notImplemented("DotProdMultiLocal");
  }

  virtual SUNErrCode DotProdMultiAllReduce(Reals /* sum */) const
  {
    notImplemented("DotProdMultiAllReduce");
  }

  // ---------------------------------------------------------------------------
  // XBraid interface operations
  // ---------------------------------------------------------------------------

  virtual SUNErrCode BufSize(sunindextype* /* size */) const
  {
    notImplemented("BufSize");
  }

  virtual SUNErrCode BufPack(void* /* buf */) const
  {
    notImplemented("BufPack");
  }

  virtual SUNErrCode BufUnpack(void* /* buf */)
  {
    notImplemented("BufUnpack");
  }

  // ---------------------------------------------------------------------------
  // Debugging output operations
  // ---------------------------------------------------------------------------

  virtual void Print() const
  {
    notImplemented("Print");
  }

  virtual void PrintFile(FILE* /* outfile */) const
  {
    notImplemented("PrintFile");
  }

  // Throws std::logic_error: the operation is not implemented by Derived. A
  // Derived member that calls an operation it does not override reaches this.
  [[noreturn]] static void notImplemented(const char* op)
  {
    throw std::logic_error(std::string("sundials::NVectorBase: operation '") +
                           op + "' is not implemented");
  }

private:
  // Converts a count passed by the C API to a view size
  static std::size_t toSize(int count) noexcept
  {
    return static_cast<std::size_t>(count);
  }

  // Vector with handle vec, as the base class, for calling operations
  static NVectorBase* base(N_Vector vec) noexcept { return FromNVector(vec); }

  static const NVectorBase* cbase(N_Vector vec) noexcept
  {
    return FromNVectorConst(vec);
  }

  // ===========================================================================
  // Static trampoline functions
  //
  // These are the glue between SUNDIALS C code and the derived C++ class. Each
  // trampoline is a static noexcept function with the signature SUNDIALS
  // expects. Trampoline nv<Method> calls the virtual function <Method> through
  // the base class, inside a guard that reports exceptions (see
  // sundials_object_base.hpp).
  // ===========================================================================

  // Category 2
  static N_Vector_ID nvGetVectorID(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetVectorID", __FILE__, __LINE__, SUNDIALS_NVEC_CUSTOM,
                            [&] { return cbase(vec)->GetVectorID(); });
  }

  // Allocates a new Derived object by calling Clone() on the source and marks
  // it as owned by its handle. The returned N_Vector is owned by SUNDIALS and
  // is freed by N_VDestroy -> nvDestroy. Every NULL return reports an error,
  // because N_VCloneVectorArray detects failure only through the last error.
  static N_Vector nvClone(N_Vector w) noexcept
  {
    return impl::guardValue(w->sunctx, "N_VClone", __FILE__, __LINE__, N_Vector{nullptr},
                            [&]() -> N_Vector {
      Derived* clone = cbase(w)->Clone(w->sunctx);
      if (clone == nullptr)
      {
        SUNContext_HandleError(w->sunctx, SUN_ERR_OP_FAIL, __LINE__,
                               "N_VClone", __FILE__, "Clone returned NULL");
        return nullptr;
      }
      if (!static_cast<NVectorBase*>(clone)->HasNVector())
      {
        delete clone;
        SUNContext_HandleError(w->sunctx, SUN_ERR_OP_FAIL, __LINE__,
                               "N_VClone", __FILE__,
                               "Clone returned an object without an N_Vector");
        return nullptr;
      }
      return Base::objAdoptClone(clone);
    });
  }

  // Always installed. Deletes the object if it was created by nvClone; for a
  // user-owned object this does nothing (its destructor frees the handle).
  static void nvDestroy(N_Vector vec) noexcept
  {
    if (vec) { Base::objDestroy(vec); }
  }

  // Category 1
  static sunindextype nvGetLength(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetLength", __FILE__, __LINE__, sunindextype{0},
                            [&] { return cbase(vec)->GetLength(); });
  }

  // z = a*x + b*y  (dispatched on z)
  static void nvLinearSum(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                          N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VLinearSum", __FILE__, __LINE__, [&] {
      base(z)->LinearSum(a, *FromNVectorConst(x), b, *FromNVectorConst(y));
    });
  }

  // z[i] = c for all i
  static void nvConst(sunrealtype c, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VConst", __FILE__, __LINE__, [&] { base(z)->Const(c); });
  }

  // z = c * x
  static void nvScale(sunrealtype c, N_Vector x, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VScale", __FILE__, __LINE__,
                    [&] { base(z)->Scale(c, *FromNVectorConst(x)); });
  }

  // Returns x . y
  static sunrealtype nvDotProd(N_Vector x, N_Vector y) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VDotProd", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->DotProd(*FromNVectorConst(y));
    });
  }

  // Returns max_i |x_i|
  static sunrealtype nvMaxNorm(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VMaxNorm", __FILE__, __LINE__, sunrealtype{0},
                            [&] { return cbase(x)->MaxNorm(); });
  }

  // ---------------------------------------------------------------------------
  // Utility operations trampolines
  // ---------------------------------------------------------------------------

  static SUNComm nvGetCommunicator(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetCommunicator", __FILE__, __LINE__, SUNComm(SUN_COMM_NULL), [&] {
      return cbase(vec)->GetCommunicator();
    });
  }

  static sunindextype nvGetLocalLength(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetLocalLength", __FILE__, __LINE__, sunindextype{0}, [&] {
      return cbase(vec)->GetLocalLength();
    });
  }

  static sunrealtype* nvGetArrayPointer(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetArrayPointer", __FILE__, __LINE__, static_cast<sunrealtype*>(nullptr), [&] {
      return base(vec)->GetArrayPointer();
    });
  }

  static void nvSetArrayPointer(sunrealtype* data, N_Vector vec) noexcept
  {
    impl::guardVoid(vec->sunctx, "N_VSetArrayPointer", __FILE__, __LINE__, [&] {
      base(vec)->SetArrayPointer(data);
    });
  }

  static sunrealtype* nvGetDeviceArrayPointer(N_Vector vec) noexcept
  {
    return impl::guardValue(vec->sunctx, "N_VGetDeviceArrayPointer", __FILE__, __LINE__, static_cast<sunrealtype*>(nullptr), [&] {
      return base(vec)->GetDeviceArrayPointer();
    });
  }

  static void nvSetDeviceArrayPointer(sunrealtype* data, N_Vector vec) noexcept
  {
    impl::guardVoid(vec->sunctx, "N_VSetDeviceArrayPointer", __FILE__, __LINE__, [&] {
      base(vec)->SetDeviceArrayPointer(data);
    });
  }

  // ---------------------------------------------------------------------------
  // Standard operations trampolines
  // ---------------------------------------------------------------------------

  static void nvProd(N_Vector x, N_Vector y, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VProd", __FILE__, __LINE__, [&] {
      base(z)->Prod(*FromNVectorConst(x), *FromNVectorConst(y));
    });
  }

  static void nvDiv(N_Vector x, N_Vector y, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VDiv", __FILE__, __LINE__, [&] {
      base(z)->Div(*FromNVectorConst(x), *FromNVectorConst(y));
    });
  }

  static void nvAbs(N_Vector x, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VAbs", __FILE__, __LINE__, [&] {
      base(z)->Abs(*FromNVectorConst(x));
    });
  }

  static void nvInv(N_Vector x, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VInv", __FILE__, __LINE__, [&] {
      base(z)->Inv(*FromNVectorConst(x));
    });
  }

  static void nvAddConst(N_Vector x, sunrealtype b, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VAddConst", __FILE__, __LINE__, [&] {
      base(z)->AddConst(*FromNVectorConst(x), b);
    });
  }

  static void nvCompare(sunrealtype c, N_Vector x, N_Vector z) noexcept
  {
    impl::guardVoid(z->sunctx, "N_VCompare", __FILE__, __LINE__, [&] {
      base(z)->Compare(c, *FromNVectorConst(x));
    });
  }

  static sunrealtype nvWrmsNorm(N_Vector x, N_Vector w) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VWrmsNorm", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->WrmsNorm(*FromNVectorConst(w));
    });
  }

  static sunrealtype nvWrmsNormMask(N_Vector x, N_Vector w, N_Vector id) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VWrmsNormMask", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->WrmsNormMask(*FromNVectorConst(w), *FromNVectorConst(id));
    });
  }

  static sunrealtype nvMin(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VMin", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->Min();
    });
  }

  static sunrealtype nvWL2Norm(N_Vector x, N_Vector w) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VWL2Norm", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->WL2Norm(*FromNVectorConst(w));
    });
  }

  static sunrealtype nvL1Norm(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VL1Norm", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->L1Norm();
    });
  }

  static sunrealtype nvMinQuotient(N_Vector num, N_Vector denom) noexcept
  {
    return impl::guardValue(num->sunctx, "N_VMinQuotient", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(num)->MinQuotient(*FromNVectorConst(denom));
    });
  }

  static sunbooleantype nvInvTest(N_Vector x, N_Vector z) noexcept
  {
    return impl::guardBool(z->sunctx, "N_VInvTest", __FILE__, __LINE__, [&] {
      return base(z)->InvTest(*FromNVectorConst(x));
    });
  }

  // Dispatched on m (the output), whereas the C implementation dispatches on x;
  // the difference is harmless for homogeneous vector types.
  static sunbooleantype nvConstrMask(N_Vector c, N_Vector x, N_Vector m) noexcept
  {
    return impl::guardBool(m->sunctx, "N_VConstrMask", __FILE__, __LINE__, [&] {
      return base(m)->ConstrMask(*FromNVectorConst(c), *FromNVectorConst(x));
    });
  }

  // ---------------------------------------------------------------------------
  // Fused operations trampolines
  // ---------------------------------------------------------------------------

  static SUNErrCode nvLinearCombination(int nvec, sunrealtype* c, N_Vector* X, N_Vector z) noexcept
  {
    return impl::guardErr(z->sunctx, "N_VLinearCombination", __FILE__, __LINE__, [&] {
      return base(z)->LinearCombination(ConstReals(c, toSize(nvec)), ConstVectors(X, toSize(nvec)));
    });
  }

  static SUNErrCode nvScaleAddMulti(int nvec, sunrealtype* a, N_Vector x, N_Vector* Y, N_Vector* Z) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VScaleAddMulti", __FILE__, __LINE__, [&] {
      return cbase(x)->ScaleAddMulti(ConstReals(a, toSize(nvec)), ConstVectors(Y, toSize(nvec)), Vectors(Z, toSize(nvec)));
    });
  }

  static SUNErrCode nvDotProdMulti(int nvec, N_Vector x, N_Vector* Y, sunrealtype* dotprods) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VDotProdMulti", __FILE__, __LINE__, [&] {
      return cbase(x)->DotProdMulti(ConstVectors(Y, toSize(nvec)), Reals(dotprods, toSize(nvec)));
    });
  }

  // ---------------------------------------------------------------------------
  // Vector array operations trampolines
  // ---------------------------------------------------------------------------

  static SUNErrCode nvLinearSumVectorArray(int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z) noexcept
  {
    return impl::guardErr(Z[0]->sunctx, "N_VLinearSumVectorArray", __FILE__, __LINE__, [&] {
      return base(Z[0])->LinearSumVectorArray(a, ConstVectors(X, toSize(nvec)), b, ConstVectors(Y, toSize(nvec)), Vectors(Z, toSize(nvec)));
    });
  }

  static SUNErrCode nvScaleVectorArray(int nvec, sunrealtype* c, N_Vector* X, N_Vector* Z) noexcept
  {
    return impl::guardErr(Z[0]->sunctx, "N_VScaleVectorArray", __FILE__, __LINE__, [&] {
      return base(Z[0])->ScaleVectorArray(ConstReals(c, toSize(nvec)), ConstVectors(X, toSize(nvec)), Vectors(Z, toSize(nvec)));
    });
  }

  static SUNErrCode nvConstVectorArray(int nvec, sunrealtype c, N_Vector* Z) noexcept
  {
    return impl::guardErr(Z[0]->sunctx, "N_VConstVectorArray", __FILE__, __LINE__, [&] {
      return base(Z[0])->ConstVectorArray(c, Vectors(Z, toSize(nvec)));
    });
  }

  static SUNErrCode nvWrmsNormVectorArray(int nvec, N_Vector* X, N_Vector* W, sunrealtype* nrm) noexcept
  {
    return impl::guardErr(X[0]->sunctx, "N_VWrmsNormVectorArray", __FILE__, __LINE__, [&] {
      return cbase(X[0])->WrmsNormVectorArray(ConstVectors(X, toSize(nvec)), ConstVectors(W, toSize(nvec)), Reals(nrm, toSize(nvec)));
    });
  }

  // Dispatched on X[0], whereas the C implementation dispatches on id; the
  // difference is harmless for homogeneous vector types.
  static SUNErrCode nvWrmsNormMaskVectorArray(int nvec, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* nrm) noexcept
  {
    return impl::guardErr(X[0]->sunctx, "N_VWrmsNormMaskVectorArray", __FILE__, __LINE__, [&] {
      return cbase(X[0])->WrmsNormMaskVectorArray(ConstVectors(X, toSize(nvec)), ConstVectors(W, toSize(nvec)), *FromNVectorConst(id), Reals(nrm, toSize(nvec)));
    });
  }

  static SUNErrCode nvScaleAddMultiVectorArray(int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z) noexcept
  {
    return impl::guardErr(X[0]->sunctx, "N_VScaleAddMultiVectorArray", __FILE__, __LINE__, [&] {
      return cbase(X[0])->ScaleAddMultiVectorArray(ConstReals(a, toSize(nsum)), ConstVectors(X, toSize(nvec)), ConstVectors2D(Y, toSize(nsum), toSize(nvec)), Vectors2D(Z, toSize(nsum), toSize(nvec)));
    });
  }

  static SUNErrCode nvLinearCombinationVectorArray(int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z) noexcept
  {
    return impl::guardErr(Z[0]->sunctx, "N_VLinearCombinationVectorArray", __FILE__, __LINE__, [&] {
      return base(Z[0])->LinearCombinationVectorArray(ConstReals(c, toSize(nsum)), ConstVectors2D(X, toSize(nsum), toSize(nvec)), Vectors(Z, toSize(nvec)));
    });
  }

  // ---------------------------------------------------------------------------
  // Local reduction operations trampolines
  // ---------------------------------------------------------------------------

  static sunrealtype nvDotProdLocal(N_Vector x, N_Vector y) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VDotProdLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->DotProdLocal(*FromNVectorConst(y));
    });
  }

  static sunrealtype nvMaxNormLocal(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VMaxNormLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->MaxNormLocal();
    });
  }

  static sunrealtype nvMinLocal(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VMinLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->MinLocal();
    });
  }

  static sunrealtype nvL1NormLocal(N_Vector x) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VL1NormLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->L1NormLocal();
    });
  }

  static sunbooleantype nvInvTestLocal(N_Vector x, N_Vector z) noexcept
  {
    return impl::guardBool(z->sunctx, "N_VInvTestLocal", __FILE__, __LINE__, [&] {
      return base(z)->InvTestLocal(*FromNVectorConst(x));
    });
  }

  // Dispatched on m (the output), whereas the C implementation dispatches on x;
  // the difference is harmless for homogeneous vector types.
  static sunbooleantype nvConstrMaskLocal(N_Vector c, N_Vector x, N_Vector m) noexcept
  {
    return impl::guardBool(m->sunctx, "N_VConstrMaskLocal", __FILE__, __LINE__, [&] {
      return base(m)->ConstrMaskLocal(*FromNVectorConst(c), *FromNVectorConst(x));
    });
  }

  static sunrealtype nvMinQuotientLocal(N_Vector num, N_Vector denom) noexcept
  {
    return impl::guardValue(num->sunctx, "N_VMinQuotientLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(num)->MinQuotientLocal(*FromNVectorConst(denom));
    });
  }

  static sunrealtype nvWSqrSumLocal(N_Vector x, N_Vector w) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VWSqrSumLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->WSqrSumLocal(*FromNVectorConst(w));
    });
  }

  static sunrealtype nvWSqrSumMaskLocal(N_Vector x, N_Vector w, N_Vector id) noexcept
  {
    return impl::guardValue(x->sunctx, "N_VWSqrSumMaskLocal", __FILE__, __LINE__, sunrealtype{0}, [&] {
      return cbase(x)->WSqrSumMaskLocal(*FromNVectorConst(w), *FromNVectorConst(id));
    });
  }

  // ---------------------------------------------------------------------------
  // Single buffer reduction operations trampolines
  // ---------------------------------------------------------------------------

  static SUNErrCode nvDotProdMultiLocal(int nvec, N_Vector x, N_Vector* Y, sunrealtype* dotprods) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VDotProdMultiLocal", __FILE__, __LINE__, [&] {
      return cbase(x)->DotProdMultiLocal(ConstVectors(Y, toSize(nvec)), Reals(dotprods, toSize(nvec)));
    });
  }

  static SUNErrCode nvDotProdMultiAllReduce(int nvec_total, N_Vector x, sunrealtype* sum) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VDotProdMultiAllReduce", __FILE__, __LINE__, [&] {
      return cbase(x)->DotProdMultiAllReduce(Reals(sum, toSize(nvec_total)));
    });
  }

  // ---------------------------------------------------------------------------
  // XBraid interface operations trampolines
  // ---------------------------------------------------------------------------

  static SUNErrCode nvBufSize(N_Vector x, sunindextype* size) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VBufSize", __FILE__, __LINE__, [&] {
      return cbase(x)->BufSize(size);
    });
  }

  static SUNErrCode nvBufPack(N_Vector x, void* buf) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VBufPack", __FILE__, __LINE__, [&] {
      return cbase(x)->BufPack(buf);
    });
  }

  static SUNErrCode nvBufUnpack(N_Vector x, void* buf) noexcept
  {
    return impl::guardErr(x->sunctx, "N_VBufUnpack", __FILE__, __LINE__, [&] {
      return base(x)->BufUnpack(buf);
    });
  }

  // ---------------------------------------------------------------------------
  // Debugging output operations trampolines
  // ---------------------------------------------------------------------------

  static void nvPrint(N_Vector vec) noexcept
  {
    impl::guardVoid(vec->sunctx, "N_VPrint", __FILE__, __LINE__, [&] { cbase(vec)->Print(); });
  }

  static void nvPrintFile(N_Vector vec, FILE* outfile) noexcept
  {
    impl::guardVoid(vec->sunctx, "N_VPrintFile", __FILE__, __LINE__, [&] { cbase(vec)->PrintFile(outfile); });
  }

  // ===========================================================================
  // Installation of the ops table
  // ===========================================================================

  // Whether to install an optional operation: true if Derived redeclares it.
  // The redeclaration must match the base declaration, so a near miss (e.g., a
  // different const or parameter type) is an error rather than a silent skip.
  template<class DerivedFn, class BaseFn>
  static constexpr bool installIfOverridden()
  {
    if constexpr (!impl::is_overridden_v<DerivedFn, BaseFn>) { return false; }
    else
    {
      static_assert(impl::same_signature_v<DerivedFn, BaseFn>,
                    "sundials::NVectorBase: Derived redeclares an operation "
                    "with a different signature than the base class "
                    "declaration; use 'override'");
      return true;
    }
  }

  // The hook called by impl::ObjectBase on every new handle.
  static void sunInstallOps(N_Vector vec) noexcept
  {
    constexpr bool derives_publicly =
      std::is_base_of_v<NVectorBase, Derived> &&
      std::is_convertible_v<Derived*, NVectorBase*>;

    static_assert(derives_publicly,
                  "sundials::NVectorBase: Derived must inherit publicly from "
                  "NVectorBase<Derived> and pass itself as the template "
                  "argument");

    // Stop after the check above if it fails, so that a bad inheritance gives
    // one clear error instead of a cascade of access errors below.
    if constexpr (derives_publicly)
    {
      N_Vector_Ops ops = vec->ops;

      // Always installed: destroy, categories 1 and 2
      ops->nvdestroy     = &nvDestroy;
      ops->nvgetvectorid = &nvGetVectorID;
      ops->nvclone       = &nvClone;
      ops->nvgetlength   = &nvGetLength;
      ops->nvlinearsum   = &nvLinearSum;
      ops->nvconst       = &nvConst;
      ops->nvscale       = &nvScale;
      ops->nvdotprod     = &nvDotProd;
      ops->nvmaxnorm     = &nvMaxNorm;

      // Category 3: installed only if Derived overrides the operation. A compile
      // error on one of the blocks below usually means one of:
      //  - Derived has another member with the op's name (an overloaded helper
      //    such as Print(std::ostream&), or a template): rename the helper;
      //  - Derived's override is private or protected: make it public, or add
      //    'friend sundials::NVectorBase<MyVec>;';
      //  - the signature differs from the base declaration: copy it and add
      //    'override'.

      // Utility operations
      if constexpr (installIfOverridden<decltype(&Derived::GetCommunicator),
                                        decltype(&NVectorBase::GetCommunicator)>())
      {
        ops->nvgetcommunicator = &nvGetCommunicator;
      }
      if constexpr (installIfOverridden<decltype(&Derived::GetLocalLength),
                                        decltype(&NVectorBase::GetLocalLength)>())
      {
        ops->nvgetlocallength = &nvGetLocalLength;
      }
      if constexpr (installIfOverridden<decltype(&Derived::GetArrayPointer),
                                        decltype(&NVectorBase::GetArrayPointer)>())
      {
        ops->nvgetarraypointer = &nvGetArrayPointer;
      }
      if constexpr (installIfOverridden<decltype(&Derived::SetArrayPointer),
                                        decltype(&NVectorBase::SetArrayPointer)>())
      {
        ops->nvsetarraypointer = &nvSetArrayPointer;
      }
      if constexpr (installIfOverridden<decltype(&Derived::GetDeviceArrayPointer),
                                        decltype(&NVectorBase::GetDeviceArrayPointer)>())
      {
        ops->nvgetdevicearraypointer = &nvGetDeviceArrayPointer;
      }
      if constexpr (installIfOverridden<decltype(&Derived::SetDeviceArrayPointer),
                                        decltype(&NVectorBase::SetDeviceArrayPointer)>())
      {
        ops->nvsetdevicearraypointer = &nvSetDeviceArrayPointer;
      }

      // Standard operations
      if constexpr (installIfOverridden<decltype(&Derived::Prod),
                                        decltype(&NVectorBase::Prod)>())
      {
        ops->nvprod = &nvProd;
      }
      if constexpr (installIfOverridden<decltype(&Derived::Div),
                                        decltype(&NVectorBase::Div)>())
      {
        ops->nvdiv = &nvDiv;
      }
      if constexpr (installIfOverridden<decltype(&Derived::Abs),
                                        decltype(&NVectorBase::Abs)>())
      {
        ops->nvabs = &nvAbs;
      }
      if constexpr (installIfOverridden<decltype(&Derived::Inv),
                                        decltype(&NVectorBase::Inv)>())
      {
        ops->nvinv = &nvInv;
      }
      if constexpr (installIfOverridden<decltype(&Derived::AddConst),
                                        decltype(&NVectorBase::AddConst)>())
      {
        ops->nvaddconst = &nvAddConst;
      }
      if constexpr (installIfOverridden<decltype(&Derived::Compare),
                                        decltype(&NVectorBase::Compare)>())
      {
        ops->nvcompare = &nvCompare;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WrmsNorm),
                                        decltype(&NVectorBase::WrmsNorm)>())
      {
        ops->nvwrmsnorm = &nvWrmsNorm;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WrmsNormMask),
                                        decltype(&NVectorBase::WrmsNormMask)>())
      {
        ops->nvwrmsnormmask = &nvWrmsNormMask;
      }
      if constexpr (installIfOverridden<decltype(&Derived::Min),
                                        decltype(&NVectorBase::Min)>())
      {
        ops->nvmin = &nvMin;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WL2Norm),
                                        decltype(&NVectorBase::WL2Norm)>())
      {
        ops->nvwl2norm = &nvWL2Norm;
      }
      if constexpr (installIfOverridden<decltype(&Derived::L1Norm),
                                        decltype(&NVectorBase::L1Norm)>())
      {
        ops->nvl1norm = &nvL1Norm;
      }
      if constexpr (installIfOverridden<decltype(&Derived::MinQuotient),
                                        decltype(&NVectorBase::MinQuotient)>())
      {
        ops->nvminquotient = &nvMinQuotient;
      }
      if constexpr (installIfOverridden<decltype(&Derived::InvTest),
                                        decltype(&NVectorBase::InvTest)>())
      {
        ops->nvinvtest = &nvInvTest;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ConstrMask),
                                        decltype(&NVectorBase::ConstrMask)>())
      {
        ops->nvconstrmask = &nvConstrMask;
      }

      // Fused operations
      if constexpr (installIfOverridden<decltype(&Derived::LinearCombination),
                                        decltype(&NVectorBase::LinearCombination)>())
      {
        ops->nvlinearcombination = &nvLinearCombination;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ScaleAddMulti),
                                        decltype(&NVectorBase::ScaleAddMulti)>())
      {
        ops->nvscaleaddmulti = &nvScaleAddMulti;
      }
      if constexpr (installIfOverridden<decltype(&Derived::DotProdMulti),
                                        decltype(&NVectorBase::DotProdMulti)>())
      {
        ops->nvdotprodmulti = &nvDotProdMulti;
      }

      // Vector array operations
      if constexpr (installIfOverridden<decltype(&Derived::LinearSumVectorArray),
                                        decltype(&NVectorBase::LinearSumVectorArray)>())
      {
        ops->nvlinearsumvectorarray = &nvLinearSumVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ScaleVectorArray),
                                        decltype(&NVectorBase::ScaleVectorArray)>())
      {
        ops->nvscalevectorarray = &nvScaleVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ConstVectorArray),
                                        decltype(&NVectorBase::ConstVectorArray)>())
      {
        ops->nvconstvectorarray = &nvConstVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WrmsNormVectorArray),
                                        decltype(&NVectorBase::WrmsNormVectorArray)>())
      {
        ops->nvwrmsnormvectorarray = &nvWrmsNormVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WrmsNormMaskVectorArray),
                                        decltype(&NVectorBase::WrmsNormMaskVectorArray)>())
      {
        ops->nvwrmsnormmaskvectorarray = &nvWrmsNormMaskVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ScaleAddMultiVectorArray),
                                        decltype(&NVectorBase::ScaleAddMultiVectorArray)>())
      {
        ops->nvscaleaddmultivectorarray = &nvScaleAddMultiVectorArray;
      }
      if constexpr (installIfOverridden<decltype(&Derived::LinearCombinationVectorArray),
                                        decltype(&NVectorBase::LinearCombinationVectorArray)>())
      {
        ops->nvlinearcombinationvectorarray = &nvLinearCombinationVectorArray;
      }

      // Local reduction operations
      if constexpr (installIfOverridden<decltype(&Derived::DotProdLocal),
                                        decltype(&NVectorBase::DotProdLocal)>())
      {
        ops->nvdotprodlocal = &nvDotProdLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::MaxNormLocal),
                                        decltype(&NVectorBase::MaxNormLocal)>())
      {
        ops->nvmaxnormlocal = &nvMaxNormLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::MinLocal),
                                        decltype(&NVectorBase::MinLocal)>())
      {
        ops->nvminlocal = &nvMinLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::L1NormLocal),
                                        decltype(&NVectorBase::L1NormLocal)>())
      {
        ops->nvl1normlocal = &nvL1NormLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::InvTestLocal),
                                        decltype(&NVectorBase::InvTestLocal)>())
      {
        ops->nvinvtestlocal = &nvInvTestLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::ConstrMaskLocal),
                                        decltype(&NVectorBase::ConstrMaskLocal)>())
      {
        ops->nvconstrmasklocal = &nvConstrMaskLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::MinQuotientLocal),
                                        decltype(&NVectorBase::MinQuotientLocal)>())
      {
        ops->nvminquotientlocal = &nvMinQuotientLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WSqrSumLocal),
                                        decltype(&NVectorBase::WSqrSumLocal)>())
      {
        ops->nvwsqrsumlocal = &nvWSqrSumLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::WSqrSumMaskLocal),
                                        decltype(&NVectorBase::WSqrSumMaskLocal)>())
      {
        ops->nvwsqrsummasklocal = &nvWSqrSumMaskLocal;
      }

      // Single buffer reduction operations
      if constexpr (installIfOverridden<decltype(&Derived::DotProdMultiLocal),
                                        decltype(&NVectorBase::DotProdMultiLocal)>())
      {
        ops->nvdotprodmultilocal = &nvDotProdMultiLocal;
      }
      if constexpr (installIfOverridden<decltype(&Derived::DotProdMultiAllReduce),
                                        decltype(&NVectorBase::DotProdMultiAllReduce)>())
      {
        ops->nvdotprodmultiallreduce = &nvDotProdMultiAllReduce;
      }

      // XBraid interface operations
      if constexpr (installIfOverridden<decltype(&Derived::BufSize),
                                        decltype(&NVectorBase::BufSize)>())
      {
        ops->nvbufsize = &nvBufSize;
      }
      if constexpr (installIfOverridden<decltype(&Derived::BufPack),
                                        decltype(&NVectorBase::BufPack)>())
      {
        ops->nvbufpack = &nvBufPack;
      }
      if constexpr (installIfOverridden<decltype(&Derived::BufUnpack),
                                        decltype(&NVectorBase::BufUnpack)>())
      {
        ops->nvbufunpack = &nvBufUnpack;
      }

      // Debugging output operations
      if constexpr (installIfOverridden<decltype(&Derived::Print),
                                        decltype(&NVectorBase::Print)>())
      {
        ops->nvprint = &nvPrint;
      }
      if constexpr (installIfOverridden<decltype(&Derived::PrintFile),
                                        decltype(&NVectorBase::PrintFile)>())
      {
        ops->nvprintfile = &nvPrintFile;
      }
    }
  }
};

} // namespace sundials

#undef SUNDIALS_NVECTOR_BASE_DEBUG_CHECK

#endif // SUNDIALS_SUNDIALS_NVECTOR_BASE_HPP
