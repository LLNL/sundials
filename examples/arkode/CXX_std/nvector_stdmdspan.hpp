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
 * An example NVector for the CXX_std examples (not an installed SUNDIALS
 * module). It wraps an arbitrary std::mdspan and supports any rank, extents,
 * layout, and accessor.
 *
 * ---------------------------------------------------------------------------
 * C++ TECHNIQUES USED IN THIS FILE
 * ---------------------------------------------------------------------------
 *
 * PARAMETER PACKS AND INDEX SEQUENCES (std::index_sequence)
 *
 *   Problem: we need to call span_[i0, i1, ..., iN-1] but the number of indices
 *   depends on the rank. Solution: std::make_index_sequence<N> generates the
 *   compile-time list {0, 1, ..., N-1}, a template captures it as a "parameter
 *   pack" (Is...), and idx[Is]... expands to idx[0], idx[1], ..., idx[N-1].
 *   Used in MultiFor::unpack and computeSize.
 *
 * FOLD EXPRESSIONS
 *
 *   (s.extent(Is) * ...) multiplies all expanded values together: s.extent(0) *
 *   s.extent(1) * ... * s.extent(N-1). Used to compute the total element count
 *   from the extents.
 *
 * IF CONSTEXPR
 *
 *   Evaluates at compile time. Only the matching branch is compiled. Used in
 *   forEach to select the right iteration strategy based on the layout type,
 *   and in the conditionally-available operator[].
 *
 * SFINAE (std::enable_if_t)
 *
 *   Makes the flat operator[] disappear at compile time for non-contiguous
 *   layouts.
 *
 * CONDITIONAL BASE CLASS (specialization on a bool)
 *
 *   NVectorBase installs an optional operation only if the derived class
 *   overrides it. A class template cannot override an operation for some
 *   template arguments and not others, so MdSpanVector derives from
 *   impl::MdSpanVectorBase<MdSpanVector, IsContiguous>. The specialization for
 *   contiguous layouts overrides GetArrayPointer; the primary template does not,
 *   so the op stays NULL.
 *
 * DECLTYPE WITH FUNCTION DECLARATIONS (no body)
 *
 *   The RightOrder and LeftOrder type aliases use a trick: a function is
 *   DECLARED (not defined) and decltype extracts the return type without ever
 *   calling it. This is a pure compile-time computation to build the
 *   MultiFor<DimOrder...> type with the right dimension ordering.
 *
 * ---------------------------------------------------------------------------
 * LAYOUT-AWARE ITERATION
 * ---------------------------------------------------------------------------
 *
 * All element-wise operations iterate over multi-indices using nested loops.
 * The NESTING ORDER of those loops is chosen to give stride-1 access in the
 * innermost loop, which is critical for cache performance:
 *
 *   layout_right (row-major):
 *     Outermost -> dim 0, innermost -> dim N-1.
 *     Last index varies fastest, matching the memory layout.
 *     Loop order determined at compile time.
 *
 *   layout_left (column-major):
 *     Outermost -> dim N-1, innermost -> dim 0.
 *     First index varies fastest, matching the memory layout.
 *     Loop order determined at compile time.
 *
 *   layout_stride (arbitrary strides):
 *     Dimensions sorted by increasing stride, so the dimension with the
 *     smallest stride (most contiguous) is the innermost loop.
 *     Order computed at construction from the mapping's runtime strides
 *     and cached for the lifetime of the object.
 *
 *   Unknown/custom layouts:
 *     Default to layout_right order (dim 0 outermost).
 *
 * The iteration order does NOT affect correctness - all N_Vector operations are
 * element-wise and order-independent. Only cache performance differs.
 *
 * In all cases, the lambda/callback receives indices in DIMENSION order (i0,
 * i1, ..., iN-1) regardless of which loop is outermost. This means
 * span_[is...] always accesses the correct element through the mapping.
 *
 * ---------------------------------------------------------------------------
 * CONTIGUOUS-ONLY FLAT ACCESS
 * ---------------------------------------------------------------------------
 *
 * For layouts that are both unique and exhaustive (layout_right, layout_left,
 * and equivalent custom layouts), the GetArrayPointer operation is installed
 * (through the conditional base class) and the flat operator[] is available
 * (through SFINAE). For non-contiguous layouts the SUNDIALS op pointer stays NULL
 * and operator[] does not exist.
 *
 * ---------------------------------------------------------------------------
 * MDSPAN REQUIREMENTS
 * ---------------------------------------------------------------------------
 *
 *   MdSpanType::element_type             -> must be sunrealtype (not const)
 *   MdSpanType::rank()                   -> constexpr, number of dimensions
 *   MdSpanType::layout_type              -> the layout policy
 *   MdSpanType::mapping_type             -> the layout mapping type
 *   MdSpanType::accessor_type            -> the accessor type
 *   span.data_handle()                   -> pointer to storage
 *   span.extent(r)                       -> size of dimension r
 *   span.extents()                       -> the extents object
 *   span.mapping()                       -> the mapping object
 *   span.accessor()                      -> the accessor object
 *   span[i, j, ...]                      -> element access (C++23)
 *   MdSpanType(ptr, mapping, accessor)   -> three-argument constructor
 *   mapping.required_span_size()         -> total memory span needed
 *   mapping.stride(r)                    -> stride of dimension r
 *                                          (only for layout_stride dispatch)
 *
 * The header includes <mdspan> and so requires C++23 (or a standard library
 * that provides std::mdspan); it fails with an #error otherwise.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_STDMDSPAN_HPP
#define SUNDIALS_NVECTOR_STDMDSPAN_HPP

#include <sundials/sundials_nvector_base.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#if __has_include(<mdspan>)
#include <mdspan>
#endif
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>

#ifndef __cpp_lib_mdspan
#error "nvector_stdmdspan.hpp requires std::mdspan (C++23 <mdspan>)"
#endif

namespace sundials {
namespace example {

// ============================================================================
// Iteration helpers
//
// The core problem: we need nested loops over N dimensions where N is a
// compile-time constant, and the loop nesting order varies by layout. These
// helpers generate the loops at compile time and call a user-provided function
// callable(i0, i1, ..., iN-1) for every multi-index.
//
// The indices are always passed to callable in DIMENSION order (i0 = dim 0, i1 = dim
// 1, etc.) regardless of which dimension is the outermost loop. This is
// achieved by storing indices in an array at their correct position (idx[dim] =
// i) and then unpacking the array into the callable using the
// index_sequence technique.
// ============================================================================

namespace impl {

// --------------------------------------------------------------------------
// Compile-time ordered iteration (MultiFor)
//
// DimOrder... specifies the loop nesting from outermost to innermost.
// Example: MultiFor<2, 1, 0> generates:
//   for (idx[2] = 0..ext(2))     <- outermost
//     for (idx[1] = 0..ext(1))
//       for (idx[0] = 0..ext(0)) <- innermost
//         callable(idx[0], idx[1], idx[2])  <- always dimension order
//
// The nest() method is recursive: each call handles one loop level and recurses
// for the next. "Level" is a compile-time counter that advances through the
// DimOrder list. When Level reaches Rank, the recursion ends and unpack()
// calls callable with the accumulated indices.
// --------------------------------------------------------------------------

template<std::size_t... DimOrder>
struct MultiFor
{
  static constexpr std::size_t Rank = sizeof...(DimOrder);

  template<class Extents, class Func>
  static void apply(const Extents& ext, Func&& callable)
  {
    std::array<std::size_t, Rank> idx{};
    nest<0>(ext, callable, idx);
  }

private:
  // The dimension ordering stored as a compile-time array.
  // order[0] = outermost dimension, order[Rank-1] = innermost.
  static constexpr std::size_t order[Rank] = {DimOrder...};

  // Recursive loop generator. At each level, loops over the dimension specified
  // by order[Level] and stores the index at idx[dim].
  template<std::size_t Level, class Extents, class Func>
  static void nest(const Extents& ext, Func& callable,
                   std::array<std::size_t, Rank>& idx)
  {
    constexpr std::size_t dim = order[Level]; // which dimension this level iterates
    for (idx[dim] = 0; idx[dim] < static_cast<std::size_t>(ext.extent(dim)); ++idx[dim])
    {
      if constexpr (Level + 1 < Rank)
        nest<Level + 1>(ext, callable, idx); // recurse to next loop level
      else
        unpack(callable, idx,
                 std::make_index_sequence<Rank>{}); // innermost: call it
    }
  }

  // Convert the index array into function arguments.
  // std::make_index_sequence<Rank>{} generates {0, 1, ..., Rank-1}.
  // The template captures this as Is..., then idx[Is]... expands to
  // idx[0], idx[1], ..., idx[Rank-1] - passing the indices to callable in
  // dimension order regardless of the loop nesting.
  template<class Func, std::size_t... Is>
  static void unpack(Func& callable, const std::array<std::size_t, Rank>& idx,
                     std::index_sequence<Is...>)
  { callable(idx[Is]...); }
};

// --------------------------------------------------------------------------
// Type aliases for standard dimension orderings
//
// A function is DECLARED (not defined) whose return type is the desired
// MultiFor<...> specialization. decltype extracts the return type without
// calling the function. This is a pure compile-time mechanism to transform an
// index_sequence into a MultiFor template parameter list.
//
// RightOrder<3> = MultiFor<0, 1, 2> (dim 0 outermost -> stride-1 for row-major)
// LeftOrder<3>  = MultiFor<2, 1, 0> (dim 2 outermost -> stride-1 for col-major)
// --------------------------------------------------------------------------

// layout_right order: 0, 1, ..., Rank-1
// make_index_sequence<3> gives index_sequence<0, 1, 2>,
// the function declaration captures Is = {0, 1, 2},
// the return type is MultiFor<0, 1, 2>.
template<std::size_t... Is>
MultiFor<Is...> makeRightOrder(std::index_sequence<Is...>);

template<std::size_t Rank>
using RightOrder = decltype(makeRightOrder(std::make_index_sequence<Rank>{}));

// layout_left order: Rank-1, ..., 1, 0
// (Rank - 1 - Is)... reverses the sequence: for Rank=3 and Is={0,1,2},
// this produces {2, 1, 0}.
template<std::size_t Rank, std::size_t... Is>
MultiFor<(Rank - 1 - Is)...> makeLeftOrder(std::index_sequence<Is...>);

template<std::size_t Rank>
using LeftOrder = decltype(makeLeftOrder<Rank>(std::make_index_sequence<Rank>{}));

// --------------------------------------------------------------------------
// Runtime ordered iteration (for layout_stride and similar)
//
// Same semantics as the compile-time version, but the dimension order is stored
// in a runtime array. The recursive call uses a runtime level counter rather
// than a template parameter, so the compiler cannot fully unroll.
// --------------------------------------------------------------------------

template<std::size_t Rank>
struct MultiForRuntime
{
  template<class Extents, class Func>
  static void apply(const Extents& ext, Func&& callable,
                    const std::array<std::size_t, Rank>& dim_order)
  {
    std::array<std::size_t, Rank> idx{};
    nest(ext, callable, idx, dim_order, 0);
  }

private:
  template<class Extents, class Func>
  static void nest(const Extents& ext, Func& callable,
                   std::array<std::size_t, Rank>& idx,
                   const std::array<std::size_t, Rank>& dim_order, std::size_t level)
  {
    const std::size_t dim = dim_order[level];
    for (idx[dim] = 0; idx[dim] < static_cast<std::size_t>(ext.extent(dim)); ++idx[dim])
    {
      if (level + 1 < Rank) nest(ext, callable, idx, dim_order, level + 1);
      else unpack(callable, idx, std::make_index_sequence<Rank>{});
    }
  }

  template<class Func, std::size_t... Is>
  static void unpack(Func& callable, const std::array<std::size_t, Rank>& idx,
                     std::index_sequence<Is...>)
  { callable(idx[Is]...); }
};

// --------------------------------------------------------------------------
// Compute dimension order sorted by increasing stride (smallest stride
// -> innermost loop -> best cache behavior).
// --------------------------------------------------------------------------

template<std::size_t Rank, class MappingType>
std::array<std::size_t, Rank> strideOrder(const MappingType& mapping)
{
  std::array<std::size_t, Rank> order;
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(),
            [&](std::size_t lhs, std::size_t rhs)
            { return mapping.stride(lhs) > mapping.stride(rhs); });
  return order;
  // order[0] = largest stride (outermost), order[Rank-1] = smallest (innermost)
}

} // namespace impl

// ============================================================================
// Optional GetArrayPointer
//
// NVectorBase installs an optional operation only if Derived overrides it, and
// a class template cannot override an operation for some template arguments
// and not others. This intermediate class sits between NVectorBase and
// MdSpanVector: the primary template declares nothing, so GetArrayPointer
// resolves to the protected default in NVectorBase and the op stays NULL. The
// specialization for contiguous layouts (unique and exhaustive mappings, e.g.,
// layout_right and layout_left) overrides it, so the op is installed.
// ============================================================================

namespace impl {

/// True if the flat pointer operations apply: the mapping is both unique and
/// exhaustive and the accessor is the default one, so data_handle() is a plain
/// sunrealtype* to the first of the contiguous elements.
template<class MdSpanType>
inline constexpr bool HasFlatArray =
  MdSpanType::mapping_type::is_always_exhaustive() &&
  MdSpanType::mapping_type::is_always_unique() &&
  std::is_same_v<typename MdSpanType::accessor_type,
                 std::default_accessor<sunrealtype>>;

template<class VectorType, bool IsContiguous>
class MdSpanVectorBase : public NVectorBase<VectorType>
{
public:
  using NVectorBase<VectorType>::NVectorBase;

  MdSpanVectorBase() = default;

protected:
  // Copy and move are protected (defaulted) to close the slicing path. Only
  // MdSpanVector, which derives from this class, uses them.
  MdSpanVectorBase(const MdSpanVectorBase&)                = default;
  MdSpanVectorBase(MdSpanVectorBase&&) noexcept            = default;
  MdSpanVectorBase& operator=(const MdSpanVectorBase&)     = default;
  MdSpanVectorBase& operator=(MdSpanVectorBase&&) noexcept = default;
};

template<class VectorType>
class MdSpanVectorBase<VectorType, true> : public NVectorBase<VectorType>
{
public:
  using NVectorBase<VectorType>::NVectorBase;

  MdSpanVectorBase() = default;

  /// Flat pointer to the data.
  sunrealtype* GetArrayPointer() override
  {
    return static_cast<VectorType*>(this)->span().data_handle();
  }

protected:
  MdSpanVectorBase(const MdSpanVectorBase&)                = default;
  MdSpanVectorBase(MdSpanVectorBase&&) noexcept            = default;
  MdSpanVectorBase& operator=(const MdSpanVectorBase&)     = default;
  MdSpanVectorBase& operator=(MdSpanVectorBase&&) noexcept = default;
};

} // namespace impl

// ============================================================================
// MdSpanVector
//
// The class is move-only: the unique_ptr member that owns the storage of a
// clone deletes the copy members. The implicit move members are noexcept (see
// the static_assert on each instantiation in the examples).
//
// Do not move out of an object that SUNDIALS still holds (e.g., one obtained
// with FromNVector from a clone, or one whose N_Vector was given to an
// integrator). Static extents cannot be reset, so a moved-from MdSpanVector
// keeps its span (and length) while its owned storage is gone. A moved-from
// object may only be destroyed or assigned to.
// ============================================================================

template<class MdSpanType>
class MdSpanVector final
  : public impl::MdSpanVectorBase<MdSpanVector<MdSpanType>,
                                  impl::HasFlatArray<MdSpanType>>
{
  using Base =
    impl::MdSpanVectorBase<MdSpanVector<MdSpanType>,
                           impl::HasFlatArray<MdSpanType>>;

public:
  using mdspan_type   = MdSpanType;
  using mapping_type  = typename mdspan_type::mapping_type;
  using accessor_type = typename mdspan_type::accessor_type;
  using layout_type   = typename mdspan_type::layout_type;

  static_assert(
    std::is_same_v<typename mdspan_type::element_type, sunrealtype>,
    "MdSpanVector requires an mdspan whose element_type is sunrealtype "
    "(a mutable, non-const element type)");

  static constexpr std::size_t Rank = mdspan_type::rank();

  // --------------------------------------------------------------------------
  // Construction
  // --------------------------------------------------------------------------

  /// Empty state: no N_Vector handle. Only available when the mdspan is
  /// default constructible (the constructor is deleted for static extents).
  MdSpanVector() = default;

  /// View constructor. Wraps the user's mdspan. The caller must keep the
  /// underlying storage alive for the lifetime of this object.
  MdSpanVector(mdspan_type user_span, SUNContext sunctx)
    : Base(sunctx),
      span_(user_span),
      len_(computeSize(user_span)),
      owned_(),
      stride_order_(computeStrideOrder(user_span))
  {}

  // --------------------------------------------------------------------------
  // Lifecycle
  // --------------------------------------------------------------------------

  /// Clones allocate fresh storage and construct a new mdspan over it.
  ///
  /// required_span_size() is used instead of the product of extents because
  /// non-contiguous layouts (e.g., layout_stride) may need more memory than the
  /// logical element count - there can be gaps between elements.
  ///
  /// The mapping and accessor are copied from the original so the clone has the
  /// same shape, layout, and element access behavior.
  ///
  /// std::make_unique<sunrealtype[]>(n) allocates n elements on the heap and
  /// wraps the pointer in a unique_ptr that automatically frees the memory when
  /// the MdSpanVector is destroyed.
  ///
  /// std::move(storage) transfers ownership of the allocation into the new
  /// MdSpanVector without copying - the source unique_ptr becomes null.
  MdSpanVector* Clone(SUNContext sunctx) const override
  {
    auto storage =
      std::make_unique<sunrealtype[]>(span_.mapping().required_span_size());
    mdspan_type new_span(storage.get(), span_.mapping(), span_.accessor());
    return new MdSpanVector(std::move(storage), new_span, sunctx);
  }

  /// Returns the number of logical elements (product of all extents).
  sunindextype GetLength() const override { return len_; }

  // --------------------------------------------------------------------------
  // Data access
  // --------------------------------------------------------------------------

  /// Access the underlying mdspan for multi-dimensional element access.
  mdspan_type& span() { return span_; }

  const mdspan_type& span() const { return span_; }

  /// Returns true if this object owns its backing storage (clone).
  bool owning() const { return owned_ != nullptr; }

  // -- Contiguous-only access (conditionally available) ---------------------
  //
  // GetArrayPointer is installed only for contiguous layouts, through the
  // conditional base class impl::MdSpanVectorBase (see above). The flat
  // operator[] below is not an N_Vector operation, so it uses SFINAE on the
  // mapping: it disappears at compile time unless the mapping is both unique
  // and exhaustive.

  /// Flat element access by linear index. Only available when the mapping is
  /// both unique and exhaustive (layout_right, layout_left, or equivalent).
  template<class MdSpan = mdspan_type,
           std::enable_if_t<impl::HasFlatArray<MdSpan>, int> = 0>
  sunrealtype& operator[](sunindextype i)
  {
    return span_.data_handle()[i];
  }

  template<class MdSpan = mdspan_type,
           std::enable_if_t<impl::HasFlatArray<MdSpan>, int> = 0>
  const sunrealtype& operator[](sunindextype i) const
  {
    return span_.data_handle()[i];
  }

  // --------------------------------------------------------------------------
  // Standard vector operations
  //
  // Each operation calls forEach() with a lambda that receives the
  // multi-indices as a parameter pack (auto... is). The lambda body uses
  // span_[is...] which expands to span_[i0, i1, ..., iN-1] - the correct mdspan
  // element access for any rank.
  //
  // [&] captures all local variables by reference. This gives the lambda
  // access to span_ (for writing), x.span_ (for reading), and any scalar
  // parameters (a, b, c).
  //
  // Example expansion for a rank-2 mdspan:
  //   forEach calls the lambda with (i0, i1)
  //   span_[is...] becomes span_[i0, i1]
  //   The mdspan's operator[] maps (i0, i1) to the correct memory offset
  // --------------------------------------------------------------------------

  void LinearSum(sunrealtype a, const MdSpanVector& x, sunrealtype b,
                 const MdSpanVector& y) override
  {
    forEach([&](auto... is)
            { span_[is...] = a * x.span_[is...] + b * y.span_[is...]; });
  }

  void Const(sunrealtype c) override
  {
    forEach([&](auto... is) { span_[is...] = c; });
  }

  void Prod(const MdSpanVector& x, const MdSpanVector& y) override
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] * y.span_[is...]; });
  }

  void Div(const MdSpanVector& x, const MdSpanVector& y) override
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] / y.span_[is...]; });
  }

  void Scale(sunrealtype c, const MdSpanVector& x) override
  {
    forEach([&](auto... is) { span_[is...] = c * x.span_[is...]; });
  }

  void Abs(const MdSpanVector& x) override
  {
    forEach([&](auto... is) { span_[is...] = std::abs(x.span_[is...]); });
  }

  void Inv(const MdSpanVector& x) override
  {
    forEach([&](auto... is) { span_[is...] = SUN_RCONST(1.0) / x.span_[is...]; });
  }

  void AddConst(const MdSpanVector& x, sunrealtype b) override
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] + b; });
  }

  sunrealtype DotProd(const MdSpanVector& y) const override
  {
    sunrealtype sum = 0;
    forEach([&](auto... is) { sum += span_[is...] * y.span_[is...]; });
    return sum;
  }

  sunrealtype MaxNorm() const override
  {
    sunrealtype max_abs = 0;
    forEach([&](auto... is)
            { max_abs = std::max(max_abs, std::abs(span_[is...])); });
    return max_abs;
  }

  sunrealtype WrmsNorm(const MdSpanVector& w) const override
  {
    return std::sqrt(WSqrSumLocal(w) / static_cast<sunrealtype>(len_));
  }

  sunrealtype WrmsNormMask(const MdSpanVector& w,
                           const MdSpanVector& id) const override
  {
    return std::sqrt(WSqrSumMaskLocal(w, id) / static_cast<sunrealtype>(len_));
  }

  sunrealtype Min() const override
  {
    sunrealtype min_val = std::numeric_limits<sunrealtype>::max();
    forEach([&](auto... is) { min_val = std::min(min_val, span_[is...]); });
    return min_val;
  }

  sunrealtype WL2Norm(const MdSpanVector& w) const override
  {
    return std::sqrt(WSqrSumLocal(w));
  }

  sunrealtype L1Norm() const override
  {
    sunrealtype sum = 0;
    forEach([&](auto... is) { sum += std::abs(span_[is...]); });
    return sum;
  }

  void Compare(sunrealtype c, const MdSpanVector& x) override
  {
    forEach([&](auto... is)
            { span_[is...] = std::abs(x.span_[is...]) >= c ? SUN_RCONST(1.0) : SUN_RCONST(0.0); });
  }

  sunbooleantype InvTest(const MdSpanVector& x) override
  {
    sunbooleantype result = SUNTRUE;
    forEach(
      [&](auto... is)
      {
        if (result == SUNFALSE) return;
        if (x.span_[is...] == SUN_RCONST(0.0))
        {
          result = SUNFALSE;
          return;
        }
        span_[is...] = SUN_RCONST(1.0) / x.span_[is...];
      });
    return result;
  }

  sunbooleantype ConstrMask(const MdSpanVector& c, const MdSpanVector& x) override
  {
    sunrealtype sum = 0;
    forEach(
      [&](auto... is)
      {
        bool test    = (std::abs(c.span_[is...]) > SUN_RCONST(1.5) &&
                        c.span_[is...] * x.span_[is...] <= SUN_RCONST(0.0)) ||
                       (std::abs(c.span_[is...]) > SUN_RCONST(0.5) &&
                        c.span_[is...] * x.span_[is...] < SUN_RCONST(0.0));
        span_[is...] = test ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
        sum += span_[is...];
      });
    return (sum < SUN_RCONST(0.5));
  }

  sunrealtype MinQuotient(const MdSpanVector& denom) const override
  {
    sunrealtype min_quot = std::numeric_limits<sunrealtype>::max();
    forEach(
      [&](auto... is)
      {
        if (denom.span_[is...] != SUN_RCONST(0.0))
          min_quot = std::min(min_quot, span_[is...] / denom.span_[is...]);
      });
    return min_quot;
  }

  // --------------------------------------------------------------------------
  // Local reduction kernels
  // --------------------------------------------------------------------------

  sunrealtype WSqrSumLocal(const MdSpanVector& w) const override
  {
    sunrealtype sum = 0;
    forEach(
      [&](auto... is)
      {
        const sunrealtype prod = span_[is...] * w.span_[is...];
        sum += prod * prod;
      });
    return sum;
  }

  sunrealtype WSqrSumMaskLocal(const MdSpanVector& w,
                               const MdSpanVector& id) const override
  {
    sunrealtype sum = 0;
    forEach(
      [&](auto... is)
      {
        if (id.span_[is...] > SUN_RCONST(0.0))
        {
          const sunrealtype prod = span_[is...] * w.span_[is...];
          sum += prod * prod;
        }
      });
    return sum;
  }

  // --------------------------------------------------------------------------
  // XBraid interface
  //
  // Pack/unpack serialize logical elements into a dense buffer. The buffer
  // size is the number of logical elements, NOT required_span_size(). The
  // iteration order ensures a deterministic packing regardless of layout.
  // --------------------------------------------------------------------------

  SUNErrCode BufSize(sunindextype* size) const override
  {
    *size = len_ * static_cast<sunindextype>(sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode BufPack(void* buf) const override
  {
    sunrealtype* buffer = static_cast<sunrealtype*>(buf);
    sunindextype offset = 0;
    // The lambda captures offset by reference and increments it for each
    // element visited. forEach iterates in a deterministic order (determined by
    // the layout), so pack/unpack visit elements in the same sequence - the
    // buffer is always consistent.
    forEach([&](auto... is) { buffer[offset++] = span_[is...]; });
    return SUN_SUCCESS;
  }

  SUNErrCode BufUnpack(void* buf) override
  {
    const sunrealtype* buffer = static_cast<const sunrealtype*>(buf);
    sunindextype offset       = 0;
    forEach([&](auto... is) { span_[is...] = buffer[offset++]; });
    return SUN_SUCCESS;
  }

  // --------------------------------------------------------------------------
  // Debug
  // --------------------------------------------------------------------------

  void Print() const override { PrintFile(stdout); }

  void PrintFile(FILE* outfile) const override
  {
    sunindextype index = 0;
    forEach(
      [&](auto... is)
      {
        std::fprintf(outfile, "  data[%ld] = " SUN_FORMAT_G "\n",
                     static_cast<long>(index++), span_[is...]);
      });
  }

private:
  mdspan_type span_{}; ///< The mdspan (view or over owned storage).
  sunindextype len_ = 0; ///< Product of all extents (logical element count).
  std::unique_ptr<sunrealtype[]> owned_; ///< Backing storage for clones (null for views).

  /// Cached dimension iteration order for layout_stride. For layout_right and
  /// layout_left the iteration order is known at compile time, so this array is
  /// unused (but must exist because the member is always declared). For
  /// layout_stride, dimensions are sorted by decreasing stride at construction
  /// time so the innermost loop has the smallest stride.
  std::array<std::size_t, Rank> stride_order_{};

  // --------------------------------------------------------------------------
  // Private clone constructor
  // --------------------------------------------------------------------------

  MdSpanVector(std::unique_ptr<sunrealtype[]> storage, mdspan_type user_span,
               SUNContext sunctx)
    : Base(sunctx),
      span_(user_span),
      len_(computeSize(user_span)),
      owned_(std::move(storage)),
      stride_order_(computeStrideOrder(user_span))
  {}

  // --------------------------------------------------------------------------
  // forEach - layout-aware dispatch
  //
  // Selects the iteration strategy based on the layout type. The decision is
  // made at compile time via if constexpr - only the matching branch is
  // compiled, so there is no runtime overhead from the other branches.
  //
  // std::forward<Func>(callable) is "perfect forwarding" - it passes the lambda
  // to the iteration helper with the same value category (lvalue or rvalue) it
  // arrived with, avoiding unnecessary copies.
  //
  // forEach is marked const because it only reads span_.extents(). The lambda
  // itself determines whether span_ elements are read or written - it captures
  // span_ from the calling method's scope, where the const qualification of
  // "this" controls mutability.
  // --------------------------------------------------------------------------

  template<class Func>
  void forEach(Func&& callable) const
  {
    if constexpr (Rank == 0) { callable(); }
    else if constexpr (std::is_same_v<layout_type, std::layout_left>)
    {
      // Column-major: dim 0 innermost (smallest stride)
      impl::LeftOrder<Rank>::apply(span_.extents(), std::forward<Func>(callable));
    }
    else if constexpr (std::is_same_v<layout_type, std::layout_stride>)
    {
      // Arbitrary strides: use cached stride-sorted order
      impl::MultiForRuntime<Rank>::apply(span_.extents(),
                                         std::forward<Func>(callable),
                                         stride_order_);
    }
    else
    {
      // layout_right and unknown layouts: dim N-1 innermost
      impl::RightOrder<Rank>::apply(span_.extents(),
                                    std::forward<Func>(callable));
    }
  }

  // --------------------------------------------------------------------------
  // Compute stride-sorted iteration order. For layout_stride, sorts dimensions
  // by decreasing stride (so the smallest stride is the innermost loop). For
  // other layouts, returns a dummy value (unused but must compile).
  // --------------------------------------------------------------------------

  static std::array<std::size_t, Rank> computeStrideOrder(const mdspan_type& user_span)
  {
    std::array<std::size_t, Rank> order{};
    if constexpr (Rank == 0) { return order; }
    else if constexpr (std::is_same_v<layout_type, std::layout_stride>)
    {
      return impl::strideOrder<Rank>(user_span.mapping());
    }
    else
    {
      // Unused for compile-time dispatch, but must be initialized.
      std::iota(order.begin(), order.end(), std::size_t{0});
      return order;
    }
  }

  // --------------------------------------------------------------------------
  // Product of all extents -> logical element count.
  //
  // Uses the index_sequence technique to expand span.extent(0), span.extent(1),
  // ..., span.extent(Rank-1) and multiplies them with a fold expression:
  // (span.extent(Is) * ...) -> span.extent(0) * ... * span.extent(N-1)
  //
  // The two-function pattern (computeSize calls computeSizeImpl with
  // make_index_sequence) exists because you can't "open" an index_sequence
  // directly - the template deduction step in computeSizeImpl converts the type
  // index_sequence<0,1,2> into the usable pack Is = {0, 1, 2}.
  // --------------------------------------------------------------------------

  template<std::size_t... Is>
  static sunindextype computeSizeImpl(const mdspan_type& span,
                                      std::index_sequence<Is...>)
  {
    return (static_cast<sunindextype>(span.extent(Is)) * ...);
  }

  static sunindextype computeSize(const mdspan_type& user_span)
  {
    if constexpr (Rank == 0) return 1;
    else return computeSizeImpl(user_span, std::make_index_sequence<Rank>{});
  }
};

} // namespace example
} // namespace sundials

#endif // SUNDIALS_NVECTOR_STDMDSPAN_HPP
