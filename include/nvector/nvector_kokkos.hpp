/* Same-rank, wrap-only Kokkos N_Vector. */
#ifndef SUNDIALS_NVECTOR_KOKKOS_HPP
#define SUNDIALS_NVECTOR_KOKKOS_HPP

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#if defined(SUNDIALS_KOKKOS_EXHAUSTIVE_LAYOUTSTRIDE_CHECK) && !defined(NDEBUG)
#include <unordered_set>
#endif
#include <sundials/sundials_nvector_base.hpp>

namespace sundials::kokkos {

template<class V>
struct view_traits
{
  using input = std::remove_cv_t<std::remove_reference_t<V>>;
  // Add Unmanaged without discarding other application-selected traits (for
  // example RandomAccess).  The Kokkos trait also handles views whose traits
  // are not exactly the default MemoryTraits<> specialization.
  using canonical_traits =
    Kokkos::MemoryTraits<input::memory_traits::flags | Kokkos::Unmanaged>;
  using managed_traits =
    Kokkos::MemoryTraits<input::memory_traits::flags & ~Kokkos::Unmanaged>;
  using view_type =
    Kokkos::View<typename input::data_type, typename input::array_layout,
                 typename input::device_type, canonical_traits>;
  using managed_type =
    Kokkos::View<typename input::data_type, typename input::array_layout,
                 typename input::device_type, managed_traits>;
  static constexpr bool is_unmanaged =
    (input::memory_traits::flags & Kokkos::Unmanaged) != 0;
};
template<class V>
using canonical_view_t = typename view_traits<V>::view_type;
template<class V>
using managed_view_t = typename view_traits<V>::managed_type;
template<class V>
inline constexpr bool is_unmanaged_view_v = view_traits<V>::is_unmanaged;

template<class V, class ExecutionSpace = typename V::execution_space>
class Vector final : public sundials::NVectorBase<Vector<V, ExecutionSpace>>
{
  using self                = Vector;
  using base                = sundials::NVectorBase<self>;
  using view_t              = V;
  using managed_type        = managed_view_t<V>;
  static constexpr int rank = view_t::rank;
  static_assert(rank > 0, "Kokkos N_Vector does not support rank zero");
  static_assert(std::is_same_v<typename view_t::non_const_value_type, sunrealtype>,
                "Kokkos N_Vector requires writable sunrealtype views");
  static_assert(is_unmanaged_view_v<view_t>,
                "Vector must use its canonical unmanaged view");

  enum class storage_mode
  {
    borrowed,
    retained
  };
  // The owner must precede the operational alias: view_ may be constructed
  // from owner_ and must never outlive that handle during destruction.
  managed_type owner_{};
  view_t view_{};
  storage_mode storage_{storage_mode::borrowed};
  ExecutionSpace exec_{};
  sunindextype length_{0};
  static constexpr bool supported_rank = rank >= 1 && rank <= 6;

  static void validate_rank()
  {
    if constexpr (rank > 6)
      throw std::invalid_argument(
        "NVECTOR_KOKKOS supports ranks 1--6; Kokkos MDRangePolicy supports "
        "only ranks below 7");
  }

  template<class W>
  static sunindextype length(const W& w)
  {
    sunindextype n = 1;
    for (int d = 0; d < rank; ++d)
    {
      if (w.extent(d) >
          static_cast<std::size_t>(std::numeric_limits<sunindextype>::max()) /
            static_cast<std::size_t>(n))
        throw std::overflow_error("Kokkos N_Vector length overflow");
      n *= static_cast<sunindextype>(w.extent(d));
    }
    return n;
  }

  template<class W>
  static void validate(const W& w, const ExecutionSpace&)
  {
    validate_rank();
    using input = std::remove_cv_t<W>;
    static_assert(input::rank == rank, "Kokkos N_Vector view rank mismatch");
    static_assert(std::is_same_v<typename input::non_const_value_type, sunrealtype>,
                  "Kokkos N_Vector requires sunrealtype values");
    static_assert(Kokkos::SpaceAccessibility<
                    ExecutionSpace, typename input::memory_space>::accessible,
                  "Kokkos N_Vector execution space cannot access the view");
    static_assert(std::is_same_v<typename view_t::array_layout, Kokkos::LayoutLeft> ||
                    std::is_same_v<typename view_t::array_layout, Kokkos::LayoutRight> ||
                    std::is_same_v<typename view_t::array_layout, Kokkos::LayoutStride>,
                  "unsupported Kokkos N_Vector layout");
    for (int d = 0; d < rank; ++d)
      if (w.extent(d) == 0)
        throw std::invalid_argument(
          "Kokkos N_Vector does not support zero extents");

    // Check both the logical length and the largest byte/element offset before
    // handing metadata to Kokkos.  This also catches LayoutStride maps whose
    // span cannot be represented by the view's size type.
    (void)length(w);
    std::size_t span = 1;
    if constexpr (std::is_same_v<typename view_t::array_layout, Kokkos::LayoutStride>)
    {
      for (int d = 0; d < rank; ++d)
      {
        const std::size_t stride = w.stride(d);
        if (stride == 0)
          throw std::invalid_argument(
            "Kokkos N_Vector rejects zero-stride mappings");
        const std::size_t n = static_cast<std::size_t>(w.extent(d) - 1);
        if (n != 0 &&
            (span > std::numeric_limits<std::size_t>::max() ||
             stride > (std::numeric_limits<std::size_t>::max() - span) / n))
          throw std::overflow_error("Kokkos N_Vector stride span overflow");
        span += n * stride;
      }
      if (span >
          static_cast<std::size_t>(std::numeric_limits<sunindextype>::max()))
        throw std::overflow_error("Kokkos N_Vector stride span overflow");

// Exhaustive checking is opt-in and debug-only because its host-side
// bookkeeping is not suitable for production construction. The normal path still
// checks a bounded prefix of the map without allocating; larger maps
// retain the documented injectivity precondition.
#if defined(SUNDIALS_KOKKOS_EXHAUSTIVE_LAYOUTSTRIDE_CHECK) && !defined(NDEBUG)
      constexpr std::size_t exhaustive_limit = 1u << 20;
      const sunindextype nlogical            = length(w);
      if (nlogical <= static_cast<sunindextype>(exhaustive_limit))
      {
        std::unordered_set<std::size_t> offsets;
        offsets.reserve(static_cast<std::size_t>(nlogical));
        for (sunindextype linear = 0; linear < nlogical; ++linear)
        {
          std::size_t q = static_cast<std::size_t>(linear), off = 0;
          for (int d = rank - 1; d >= 0; --d)
          {
            const std::size_t i = q % static_cast<std::size_t>(w.extent(d));
            q /= static_cast<std::size_t>(w.extent(d));
            off += i * static_cast<std::size_t>(w.stride(d));
          }
          if (!offsets.insert(off).second)
            throw std::invalid_argument(
              "Kokkos N_Vector rejects overlapping LayoutStride mapping");
        }
      }
      else
#endif
      {
        constexpr std::size_t structural_limit = 4096;
        const sunindextype nlogical            = length(w);
        const std::size_t checked =
          nlogical <= static_cast<sunindextype>(structural_limit)
            ? static_cast<std::size_t>(nlogical)
            : structural_limit;
        for (std::size_t linear = 0; linear < checked; ++linear)
        {
          std::size_t q   = linear;
          std::size_t off = 0;
          for (int d = rank - 1; d >= 0; --d)
          {
            const std::size_t i = q % static_cast<std::size_t>(w.extent(d));
            q /= static_cast<std::size_t>(w.extent(d));
            off += i * static_cast<std::size_t>(w.stride(d));
          }
          for (std::size_t prior = 0; prior < linear; ++prior)
          {
            std::size_t p        = prior;
            std::size_t previous = 0;
            for (int d = rank - 1; d >= 0; --d)
            {
              const std::size_t i = p % static_cast<std::size_t>(w.extent(d));
              p /= static_cast<std::size_t>(w.extent(d));
              previous += i * static_cast<std::size_t>(w.stride(d));
            }
            if (off == previous)
              throw std::invalid_argument(
                "Kokkos N_Vector rejects overlapping LayoutStride mapping");
          }
        }
      }
    }
    (void)span;
  }

  void check(const self& x) const
  {
    for (int d = 0; d < rank; ++d)
      if (view_.extent(d) != x.view_.extent(d))
        throw std::invalid_argument("Kokkos N_Vector shape mismatch");
  }

  template<class F>
  void each(F f) const
  {
    if constexpr (rank == 1)
    {
      auto a = view_;
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i) { f(a(i)); });
    }
    else if constexpr (rank <= 6)
    {
      Kokkos::Array<sunindextype, rank> lo{}, hi{};
      for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
      auto a = view_;
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<rank>>(exec_, lo, hi),
        KOKKOS_LAMBDA(auto... i) { f(a(i...)); });
    }
  }

  template<class F>
  void binary(const self& x, F f)
  {
    check(x);
    auto z = view_;
    auto a = x.view_;
    if constexpr (rank == 1)
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i) { f(z(i), a(i)); });
    else if constexpr (rank <= 6)
    {
      Kokkos::Array<sunindextype, rank> lo{}, hi{};
      for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<rank>>(exec_, lo, hi),
        KOKKOS_LAMBDA(auto... i) { f(z(i...), a(i...)); });
    }
  }
  // Kokkos supports variadic MDRange indices.  These reducers use
  // rank-specific functors because the accumulator is a trailing parameter,
  // making generic parameter-pack deduction problematic.
  template<int R, class VV, class F, class T>
  struct md_unary_reducer;
#define SUNDIALS_KOKKOS_UNARY_REDUCER(R, ARGS, INDEX)             \
  template<class VV, class F, class T>                            \
  struct md_unary_reducer<R, VV, F, T>                            \
  {                                                               \
    VV a;                                                         \
    F f;                                                          \
    KOKKOS_FUNCTION void operator() ARGS const { f(a INDEX, v); } \
  };
  SUNDIALS_KOKKOS_UNARY_REDUCER(2, (sunindextype i, sunindextype j, T& v), (i, j))
  SUNDIALS_KOKKOS_UNARY_REDUCER(
    3, (sunindextype i, sunindextype j, sunindextype k, T& v), (i, j, k))
  SUNDIALS_KOKKOS_UNARY_REDUCER(4,
                                (sunindextype i, sunindextype j, sunindextype k,
                                 sunindextype l, T& v),
                                (i, j, k, l))
  SUNDIALS_KOKKOS_UNARY_REDUCER(5,
                                (sunindextype i, sunindextype j, sunindextype k,
                                 sunindextype l, sunindextype m, T& v),
                                (i, j, k, l, m))
  SUNDIALS_KOKKOS_UNARY_REDUCER(6,
                                (sunindextype i, sunindextype j, sunindextype k,
                                 sunindextype l, sunindextype m, sunindextype n,
                                 T& v),
                                (i, j, k, l, m, n))
#undef SUNDIALS_KOKKOS_UNARY_REDUCER

  template<class VV, class F, class T>
  sunrealtype md_reduce(VV a, F f, T& result) const
  {
    Kokkos::Array<sunindextype, rank> lo{}, hi{};
    for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
    Kokkos::parallel_reduce("sundials_nvector_kokkos",
                            Kokkos::MDRangePolicy<ExecutionSpace,
                                                  Kokkos::Rank<rank>>(exec_, lo,
                                                                      hi),
                            md_unary_reducer<rank, V, F, T>{a, f}, result);
    return result;
  }

  template<class VV, class F, class R>
  void md_reduce_with_reducer(VV a, F f, R& result) const
  {
    Kokkos::Array<sunindextype, rank> lo{}, hi{};
    for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
    Kokkos::parallel_reduce("sundials_nvector_kokkos",
                            Kokkos::MDRangePolicy<ExecutionSpace,
                                                  Kokkos::Rank<rank>>(exec_, lo,
                                                                      hi),
                            md_unary_reducer<rank, V, F, sunrealtype>{a, f},
                            result);
  }

  template<class A, class B, class F, class T>
  void md_reduce_binary(A a, B b, F f, T& result) const
  {
    Kokkos::Array<sunindextype, rank> lo{}, hi{};
    for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
    Kokkos::parallel_reduce("sundials_nvector_kokkos",
                            Kokkos::MDRangePolicy<ExecutionSpace,
                                                  Kokkos::Rank<rank>>(exec_, lo,
                                                                      hi),
                            md_binary_reducer<rank, A, B, F, T>{a, b, f}, result);
  }

  template<class A, class B, class F, class R>
  void md_reduce_binary_with_reducer(A a, B b, F f, R& result) const
  {
    Kokkos::Array<sunindextype, rank> lo{}, hi{};
    for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
    Kokkos::parallel_reduce("sundials_nvector_kokkos",
                            Kokkos::MDRangePolicy<ExecutionSpace,
                                                  Kokkos::Rank<rank>>(exec_, lo,
                                                                      hi),
                            md_binary_reducer<rank, A, B, F, sunrealtype>{a, b, f},
                            result);
  }
  template<int R, class A, class B, class C, class F, class T>
  struct md_ternary_reducer;

  template<class T>
  struct reduction_value
  {
    using type = T;
  };

  template<class T>
  struct reduction_value<Kokkos::LAnd<T>>
  {
    using type = T;
  };

#define SUNDIALS_KOKKOS_TERNARY_REDUCER(R, ARGS, INDEX) \
  template<class A, class B, class C, class F, class T> \
  struct md_ternary_reducer<R, A, B, C, F, T>           \
  {                                                     \
    A a;                                                \
    B b;                                                \
    C c;                                                \
    F f;                                                \
    KOKKOS_FUNCTION void operator() ARGS const          \
    { f(a INDEX, b INDEX, c INDEX, v); }                \
  };
  SUNDIALS_KOKKOS_TERNARY_REDUCER(2, (sunindextype i, sunindextype j, T& v),
                                  (i, j))
  SUNDIALS_KOKKOS_TERNARY_REDUCER(
    3, (sunindextype i, sunindextype j, sunindextype k, T& v), (i, j, k))
  SUNDIALS_KOKKOS_TERNARY_REDUCER(4,
                                  (sunindextype i, sunindextype j,
                                   sunindextype k, sunindextype l, T& v),
                                  (i, j, k, l))
  SUNDIALS_KOKKOS_TERNARY_REDUCER(5,
                                  (sunindextype i, sunindextype j, sunindextype k,
                                   sunindextype l, sunindextype m, T& v),
                                  (i, j, k, l, m))
  SUNDIALS_KOKKOS_TERNARY_REDUCER(6,
                                  (sunindextype i, sunindextype j,
                                   sunindextype k, sunindextype l,
                                   sunindextype m, sunindextype n, T& v),
                                  (i, j, k, l, m, n))
#undef SUNDIALS_KOKKOS_TERNARY_REDUCER

  template<class A, class B, class C, class F, class T>
  void md_reduce_ternary(A a, B b, C c, F f, T& r) const
  {
    Kokkos::Array<sunindextype, rank> lo{}, hi{};
    for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
    using value_type = typename reduction_value<T>::type;
    Kokkos::parallel_reduce("sundials_nvector_kokkos",
                            Kokkos::MDRangePolicy<ExecutionSpace,
                                                  Kokkos::Rank<rank>>(exec_, lo,
                                                                      hi),
                            md_ternary_reducer<rank, A, B, C, F, value_type>{a, b,
                                                                             c, f},
                            r);
  }
  template<int R, class A, class B, class F, class T>
  struct md_binary_reducer;
#define SUNDIALS_KOKKOS_BINARY_REDUCER(R, ARGS, INDEX)                     \
  template<class A, class B, class F, class T>                             \
  struct md_binary_reducer<R, A, B, F, T>                                  \
  {                                                                        \
    A a;                                                                   \
    B b;                                                                   \
    F f;                                                                   \
    KOKKOS_FUNCTION void operator() ARGS const { f(a INDEX, b INDEX, v); } \
  };
  SUNDIALS_KOKKOS_BINARY_REDUCER(2, (sunindextype i, sunindextype j, T& v),
                                 (i, j))
  SUNDIALS_KOKKOS_BINARY_REDUCER(
    3, (sunindextype i, sunindextype j, sunindextype k, T& v), (i, j, k))
  SUNDIALS_KOKKOS_BINARY_REDUCER(4,
                                 (sunindextype i, sunindextype j,
                                  sunindextype k, sunindextype l, T& v),
                                 (i, j, k, l))
  SUNDIALS_KOKKOS_BINARY_REDUCER(5,
                                 (sunindextype i, sunindextype j, sunindextype k,
                                  sunindextype l, sunindextype m, T& v),
                                 (i, j, k, l, m))
  SUNDIALS_KOKKOS_BINARY_REDUCER(6,
                                 (sunindextype i, sunindextype j,
                                  sunindextype k, sunindextype l,
                                  sunindextype m, sunindextype n, T& v),
                                 (i, j, k, l, m, n))
#undef SUNDIALS_KOKKOS_BINARY_REDUCER

  template<class F>
  void ternary(const self& x, const self& y, F f)
  {
    check(x);
    check(y);
    auto z = view_;
    auto a = x.view_;
    auto b = y.view_;
    if constexpr (rank == 1)
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i) { f(z(i), a(i), b(i)); });
    else if constexpr (rank <= 6)
    {
      Kokkos::Array<sunindextype, rank> lo{}, hi{};
      for (int d = 0; d < rank; ++d) hi[d] = view_.extent(d);
      Kokkos::parallel_for(
        "sundials_nvector_kokkos",
        Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<rank>>(exec_, lo, hi),
        KOKKOS_LAMBDA(auto... i) { f(z(i...), a(i...), b(i...)); });
    }
  }

  template<class F>
  sunrealtype reduce(F f, sunrealtype init = 0) const
  {
    sunrealtype r = init;
    auto a        = view_;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "sundials_nvector_kokkos",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) { f(a(i), v); }, r);
    else if constexpr (rank <= 6)
    {
      md_reduce(a, f, r);
    }
    return r;
  }

  sunrealtype reduce_max_abs() const
  {
    sunrealtype r = -std::numeric_limits<sunrealtype>::max();
    auto a        = view_;
    Kokkos::Max<sunrealtype> red(r);
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VMaxNorm",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) {
          v = Kokkos::max(v, std::abs(a(i)));
        },
        red);
    else if constexpr (rank <= 6)
    {
      md_reduce_with_reducer(
        a, [] KOKKOS_FUNCTION(const auto& x, sunrealtype& v)
        { v = Kokkos::max(v, std::abs(x)); }, red);
    }
    return red.reference();
  }

  sunrealtype reduce_min() const
  {
    sunrealtype r = std::numeric_limits<sunrealtype>::max();
    auto a        = view_;
    Kokkos::Min<sunrealtype> red(r);
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VMin", Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) { v = Kokkos::min(v, a(i)); }, red);
    else if constexpr (rank <= 6)
    {
      md_reduce_with_reducer(
        a, [] KOKKOS_FUNCTION(const auto& x, sunrealtype& v)
        { v = Kokkos::min(v, x); }, red);
    }
    return red.reference();
  }

  sunrealtype square_sum(const self& w) const
  {
    check(w);
    auto a        = view_;
    auto b        = w.view_;
    sunrealtype r = 0;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VWSqrSum",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) {
          v += (a(i) * b(i)) * (a(i) * b(i));
        },
        r);
    else if constexpr (rank <= 6)
      md_reduce_binary(
        a, b, [] KOKKOS_FUNCTION(const auto& x, const auto& y, sunrealtype& v)
        { v += (x * y) * (x * y); }, r);
    return r;
  }

  template<class A, class B, class M>
  void for_each_reduce(A a, B b, M m, sunrealtype& r, bool mask) const
  {
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VMaskedSum",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) {
          if (!mask || m(i) > 0) v += (a(i) * b(i)) * (a(i) * b(i));
        },
        r);
    else if constexpr (rank <= 6)
    {
      md_reduce_ternary(
        a, b, m,
        [=] KOKKOS_FUNCTION(const auto& x, const auto& y, const auto& z,
                            sunrealtype& v)
        {
          if (!mask || z > 0) v += (x * y) * (x * y);
        },
        r);
    }
  }

  sunrealtype min_quotient(const self& d) const
  {
    auto a        = view_;
    auto b        = d.view_;
    sunrealtype r = SUN_BIG_REAL;
    Kokkos::Min<sunrealtype> red(r);
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VMinQuotient",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) {
          if (b(i) != 0) v = Kokkos::min(v, a(i) / b(i));
        },
        red);
    else if constexpr (rank <= 6)
      md_reduce_binary_with_reducer(
        a, b,
        [] KOKKOS_FUNCTION(const auto& x, const auto& y, sunrealtype& v)
        {
          if (y != 0) v = Kokkos::min(v, x / y);
        },
        red);
    return red.reference();
  }

  sunbooleantype inverse_test(const self& x) const
  {
    bool ok  = true;
    auto out = view_;
    auto a   = x.view_;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VInvTest",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, bool& v) {
          if (a(i) == 0) v = false;
          else out(i) = 1 / a(i);
        },
        Kokkos::LAnd<bool>(ok));
    else if constexpr (rank <= 6)
    {
      Kokkos::LAnd<bool> red(ok);
      md_reduce_ternary(
        out, a, a,
        [] KOKKOS_FUNCTION(auto& o, const auto& z, const auto&, bool& v)
        {
          if (z == 0) v = false;
          else o = 1 / z;
        },
        red);
      ok = red.reference();
    }
    return ok ? SUNTRUE : SUNFALSE;
  }

  sunbooleantype constraint_mask(const self& c, const self& x) const
  {
    bool ok  = true;
    auto out = view_;
    auto q   = c.view_;
    auto b   = x.view_;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VConstrMask",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, bool& v) {
          const auto z   = q(i);
          const auto y   = b(i);
          const bool bad = (z == 2 && y <= 0) || (z == 1 && y < 0) ||
                           (z == -1 && y > 0) || (z == -2 && y >= 0);
          out(i)         = bad ? 1 : 0;
          if (bad) v = false;
        },
        Kokkos::LAnd<bool>(ok));
    else if constexpr (rank <= 6)
    {
      Kokkos::LAnd<bool> red(ok);
      md_reduce_ternary(
        out, q, b,
        [] KOKKOS_FUNCTION(auto& o, const auto& z, const auto& y, bool& v)
        {
          const bool bad = (z == 2 && y <= 0) || (z == 1 && y < 0) ||
                           (z == -1 && y > 0) || (z == -2 && y >= 0);
          o              = bad ? 1 : 0;
          if (bad) v = false;
        },
        red);
      ok = red.reference();
    }
    return ok ? SUNTRUE : SUNFALSE;
  }

  explicit Vector(const managed_type& v, SUNContext c, ExecutionSpace e, bool own)
    : base(c),
      owner_(v),
      view_(view_t(owner_)),
      exec_(std::move(e)),
      length_(length(v)),
      storage_(own ? storage_mode::retained : storage_mode::borrowed)
  { validate_rank(); }

public:
  using view_type         = view_t;
  using managed_view_type = managed_type;
  using execution_space   = ExecutionSpace;
  using memory_space      = typename view_t::memory_space;
  using size_type         = typename view_t::size_type;

  template<class W>
  explicit Vector(const W& w, SUNContext c, ExecutionSpace e = {})
    : base(c), owner_{}, view_(w), exec_(std::move(e)), length_(0)
  {
    static_assert(!std::is_const_v<typename W::value_type>,
                  "Kokkos N_Vector requires a writable view");
    validate(w, exec_);
    length_ = length(w);
    if constexpr (!is_unmanaged_view_v<W> &&
                  std::is_constructible_v<managed_type, W>)
    {
      // Keep the original managed handle.  Never manufacture an owner from
      // view_, which is deliberately an unmanaged operational alias.
      owner_   = managed_type(w);
      view_    = view_t(owner_);
      storage_ = storage_mode::retained;
    }
    else if constexpr (!is_unmanaged_view_v<W> &&
                       !std::is_constructible_v<managed_type, W>)
    {
      throw std::invalid_argument("Kokkos N_Vector view cannot be cloned");
    }
  }

  Vector(const Vector&)                = delete;
  Vector& operator=(const Vector&)     = delete;
  Vector(Vector&&) noexcept            = default;
  Vector& operator=(Vector&&) noexcept = default;

  view_type View() const noexcept { return view_; }

  sunindextype GetLength() const override { return length_; }

  sunindextype GetLocalLength() const override { return length_; }

  SUNComm GetCommunicator() const override { return SUN_COMM_NULL; }

  N_Vector_ID GetVectorID() const override { return SUNDIALS_NVEC_KOKKOS; }

  sunrealtype* GetDeviceArrayPointer() override
  { return rank == 1 && view_.span_is_contiguous() ? view_.data() : nullptr; }

  sunrealtype* GetArrayPointer() override
  {
    using ms = typename view_t::memory_space;
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace, ms>::accessible)
      return rank == 1 && view_.span_is_contiguous() ? view_.data() : nullptr;
    return nullptr;
  }

  self* Clone(SUNContext c) const override
  {
    managed_type w;
    if constexpr (std::is_same_v<typename view_t::array_layout, Kokkos::LayoutStride>)
    {
      Kokkos::LayoutStride l;
      for (int d = 0; d < rank; ++d)
      {
        l.dimension[d] = view_.extent(d);
        l.stride[d]    = view_.stride(d);
      }
      w = managed_type("N_Vector clone", l);
    }
    else
    {
      if constexpr (rank == 1)
        w = managed_type("N_Vector clone", view_.extent(0));
      else if constexpr (rank == 2)
        w = managed_type("N_Vector clone", view_.extent(0), view_.extent(1));
      else if constexpr (rank == 3)
        w = managed_type("N_Vector clone", view_.extent(0), view_.extent(1),
                         view_.extent(2));
      else if constexpr (rank == 4)
        w = managed_type("N_Vector clone", view_.extent(0), view_.extent(1),
                         view_.extent(2), view_.extent(3));
      else if constexpr (rank == 5)
        w = managed_type("N_Vector clone", view_.extent(0), view_.extent(1),
                         view_.extent(2), view_.extent(3), view_.extent(4));
      else if constexpr (rank == 6)
        w = managed_type("N_Vector clone", view_.extent(0), view_.extent(1),
                         view_.extent(2), view_.extent(3), view_.extent(4),
                         view_.extent(5));
       else throw std::invalid_argument("unsupported Kokkos rank");
    }
    return new self(w, c, exec_, true);
  }

  void Const(sunrealtype c) override
  {
    auto z = view_;
    each([=] KOKKOS_FUNCTION(auto& v) { v = c; });
  }

  void LinearSum(sunrealtype a, const self& x, sunrealtype b, const self& y) override
  {
    ternary(x, y, [=] KOKKOS_FUNCTION(auto& z, const auto& p, const auto& q)
            { z = a * p + b * q; });
  }

  void Scale(sunrealtype c, const self& x) override
  {
    binary(x, [=] KOKKOS_FUNCTION(auto& z, const auto& a) { z = c * a; });
  }

  sunrealtype DotProd(const self& y) const override
  {
    check(y);
    auto a        = view_;
    auto b        = y.view_;
    sunrealtype r = 0;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VDotProd",
        Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) { v += a(i) * b(i); }, r);
    else if constexpr (rank <= 6)
      md_reduce_binary(
        a, b, [] KOKKOS_FUNCTION(const auto& x, const auto& y, sunrealtype& v)
        { v += x * y; }, r);
    return r;
  }

  sunrealtype MaxNorm() const override { return reduce_max_abs(); }

  void Prod(const self& x, const self& y) override
  {
    ternary(x, y, [] KOKKOS_FUNCTION(auto& z, const auto& a, const auto& b)
            { z = a * b; });
  }

  void Div(const self& x, const self& y) override
  {
    ternary(x, y, [] KOKKOS_FUNCTION(auto& z, const auto& a, const auto& b)
            { z = a / b; });
  }

  void Abs(const self& x) override
  {
    binary(x, [] KOKKOS_FUNCTION(auto& z, const auto& a) { z = std::abs(a); });
  }

  void Inv(const self& x) override
  {
    binary(x, [] KOKKOS_FUNCTION(auto& z, const auto& a) { z = 1 / a; });
  }

  void AddConst(const self& x, sunrealtype b) override
  {
    binary(x, [=] KOKKOS_FUNCTION(auto& z, const auto& a) { z = a + b; });
  }

  void Compare(sunrealtype c, const self& x) override
  {
    binary(x, [=] KOKKOS_FUNCTION(auto& z, const auto& a)
           { z = std::abs(a) >= c ? 1 : 0; });
  }

  sunrealtype Min() const override { return reduce_min(); }

  sunrealtype L1Norm() const override
  {
    return reduce([] KOKKOS_FUNCTION(const auto& a, auto& r)
                  { r += std::abs(a); });
  }

  sunrealtype WL2Norm(const self& w) const override
  {
    check(w);
    auto a        = view_;
    auto b        = w.view_;
    sunrealtype r = 0;
    if constexpr (rank == 1)
      Kokkos::parallel_reduce(
        "N_VWL2", Kokkos::RangePolicy<ExecutionSpace>(exec_, 0, view_.extent(0)),
        KOKKOS_LAMBDA(auto i, sunrealtype& v) {
          v += (a(i) * b(i)) * (a(i) * b(i));
        },
        r);
    else if constexpr (rank <= 6)
      md_reduce_binary(
        a, b, [] KOKKOS_FUNCTION(const auto& x, const auto& y, sunrealtype& v)
        { v += (x * y) * (x * y); }, r);
    return std::sqrt(r);
  }

  sunrealtype WrmsNorm(const self& w) const override
  { return std::sqrt(WSqrSumLocal(w) / static_cast<sunrealtype>(length_)); }

  sunrealtype WrmsNormMask(const self& w, const self& id) const override
  {
    check(w);
    check(id);
    auto a        = view_;
    auto b        = w.view_;
    auto m        = id.view_;
    sunrealtype r = 0;
    for_each_reduce(a, b, m, r, true);
    return std::sqrt(r / static_cast<sunrealtype>(length_));
  }

  sunrealtype MinQuotient(const self& d) const override
  {
    check(d);
    return min_quotient(d);
  }

  sunbooleantype InvTest(const self& x) override
  {
    check(x);
    return inverse_test(x);
  }

  sunbooleantype ConstrMask(const self& c, const self& x) override
  {
    check(c);
    check(x);
    return constraint_mask(c, x);
  }

  sunrealtype WSqrSumLocal(const self& w) const override
  {
    check(w);
    return square_sum(w);
  }

  sunrealtype WSqrSumMaskLocal(const self& w, const self& id) const override
  {
    check(w);
    check(id);
    sunrealtype r = 0;
    auto a        = view_;
    auto b        = w.view_;
    auto m        = id.view_;
    for_each_reduce(a, b, m, r, true);
    return r;
  }

  sunrealtype DotProdLocal(const self& y) const override { return DotProd(y); }

  sunrealtype MaxNormLocal() const override { return MaxNorm(); }

  sunrealtype MinLocal() const override { return Min(); }

  sunrealtype L1NormLocal() const override { return L1Norm(); }

  sunbooleantype InvTestLocal(const self& x) override { return InvTest(x); }

  sunbooleantype ConstrMaskLocal(const self& c, const self& x) override
  { return ConstrMask(c, x); }

  sunrealtype MinQuotientLocal(const self& d) const override
  { return MinQuotient(d); }
};

template<class V>
using VectorFor =
  Vector<canonical_view_t<V>, typename std::remove_cv_t<V>::execution_space>;
template<class V>
Vector(V, SUNContext, typename V::execution_space = {})
  -> Vector<canonical_view_t<V>, typename V::execution_space>;
} // namespace sundials::kokkos
#endif
