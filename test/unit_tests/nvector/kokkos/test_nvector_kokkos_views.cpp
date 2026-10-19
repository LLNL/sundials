#include <Kokkos_Core.hpp>
#include <cstdio>
#include <string>
#include <stdexcept>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>

using ExecSpace = Kokkos::DefaultExecutionSpace;
using V1 = Kokkos::View<sunrealtype*, typename ExecSpace::memory_space>;
using V2 = Kokkos::View<sunrealtype**, Kokkos::LayoutRight,
                        typename ExecSpace::memory_space>;
using VL = Kokkos::View<sunrealtype**, Kokkos::LayoutLeft,
                        typename ExecSpace::memory_space>;
using N1 = sundials::kokkos::VectorFor<V1>;
using N2 = sundials::kokkos::VectorFor<V2>;
using NL = sundials::kokkos::VectorFor<VL>;
using VS = Kokkos::View<sunrealtype*, Kokkos::LayoutStride,
                        typename ExecSpace::memory_space>;
using NS = sundials::kokkos::VectorFor<VS>;
using U1 = Kokkos::View<sunrealtype*, typename ExecSpace::memory_space,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using V3 = Kokkos::View<sunrealtype***, typename ExecSpace::memory_space>;
using N3 = sundials::kokkos::VectorFor<V3>;
using V6 = Kokkos::View<sunrealtype******, typename ExecSpace::memory_space>;
using N6 = sundials::kokkos::VectorFor<V6>;
using V7 = Kokkos::View<sunrealtype*******, typename ExecSpace::memory_space>;
using V8 = Kokkos::View<sunrealtype********, typename ExecSpace::memory_space>;
using N7 = sundials::kokkos::VectorFor<V7>;
using N8 = sundials::kokkos::VectorFor<V8>;

static bool check(bool condition, const char* message)
{
  if (!condition) std::fprintf(stderr, "FAIL: %s\n", message);
  return condition;
}

#define CHECK(condition) \
  do {                    \
    if (!check((condition), #condition)) return 1; \
  } while (false)

int main(int argc, char** argv)
{
  Kokkos::initialize(argc, argv);
  {
    sundials::Context ctx;
    V1 a("a", 16);
    N1 x{a, ctx};
    U1 borrowed(a.data(), a.extent(0));
    N1 borrowed_vector{borrowed, ctx};
    CHECK(borrowed_vector.View().data() == a.data());
    CHECK(N1::FromNVector(x) == &x);
    N_Vector clone_handle = N_VClone(x);
    auto clone = N1::FromNVector(clone_handle);
    CHECK(clone->View().data() != x.View().data());
    N_VDestroy(clone_handle);

    V2 b("b", 4, 4);
    N2 y{b, ctx};
    N_VConst(SUN_RCONST(2.0), y);
    CHECK(y.GetLength() == 16);
    CHECK(y.View().extent(0) == 4 && y.View().extent(1) == 4);

    VL left("left", 4, 4);
    NL left_vector{left, ctx};
    N_VConst(SUN_RCONST(3.0), left_vector);
    CHECK(left_vector.View().stride(0) == 1);

    V3 p("p", 2, 2, 2), q("q", 2, 2, 2), out("out", 2, 2, 2);
    N3 np{p, ctx}, nq{q, ctx}, nout{out, ctx};
      // MinQuotient computes out / denominator; use distinct values so the
      // numerator and denominator roles are explicit.
      Kokkos::deep_copy(p, SUN_RCONST(4.0));
      Kokkos::deep_copy(out, SUN_RCONST(2.0));
    Kokkos::deep_copy(q, SUN_RCONST(2.0));
      CHECK(nout.MinQuotient(np) == SUN_RCONST(0.5));
    Kokkos::deep_copy(q, SUN_RCONST(0.0));
    CHECK(nout.MinQuotient(nq) == SUN_BIG_REAL);
    Kokkos::deep_copy(q, SUN_RCONST(2.0));
    CHECK(nout.InvTest(nq) == SUNTRUE);
    auto oh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);
    CHECK(oh(0, 0, 0) == SUN_RCONST(0.5));
    Kokkos::deep_copy(p, SUN_RCONST(1.0));
    CHECK(nout.ConstrMask(np, nq) == SUNTRUE);
    auto maskh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);
    CHECK(maskh(0, 0, 0) == SUN_RCONST(0.0));

    // Rank six is the highest native MDRangePolicy supported by Kokkos.
    V6 rank6("rank6", 1, 1, 1, 1, 1, 2);
    N6 n6{rank6, ctx};
    N_VConst(SUN_RCONST(4.0), n6);
    auto rank6_clone_handle = N_VClone(n6);
    auto rank6_clone = N6::FromNVector(rank6_clone_handle);
    CHECK(rank6_clone->View().rank == 6);
    CHECK(rank6_clone->View().extent(5) == 2);
    CHECK(rank6_clone->View().data() != n6.View().data());
    N_VDestroy(rank6_clone_handle);

    // High ranks are rejected before ownership is retained or work is queued.
    for (int rejected_rank : {7, 8})
    {
      bool rank_rejected = false;
      try
      {
        if (rejected_rank == 7)
        {
          V7 rank7("rank7", 1, 1, 1, 1, 1, 1, 1);
          N7 invalid{rank7, ctx};
        }
        else
        {
          V8 rank8("rank8", 1, 1, 1, 1, 1, 1, 1, 1);
          // Construction must fail without launching an operation.
          N8 invalid{rank8, ctx};
        }
      }
      catch (const std::invalid_argument& error)
      {
        rank_rejected = std::string(error.what()).find("ranks 1--6") !=
                        std::string::npos;
      }
      CHECK(rank_rejected);
    }

    // A subview is an ordinary strided view; clones retain its mapping and
    // remain independent storage.
    V1 backing("backing", 32);
    auto sub = Kokkos::subview(backing, Kokkos::pair<std::size_t, std::size_t>(1, 17));
    NS s{sub, ctx};
    auto sh = N_VClone(s);
    auto sc = NS::FromNVector(sh);
    CHECK(sc->View().extent(0) == s.View().extent(0));
    CHECK(sc->View().stride(0) == s.View().stride(0));
    CHECK(sc->View().data() != s.View().data());
    N_VDestroy(sh);

    bool rejected = false;
    try {
      VS zero("zero", Kokkos::LayoutStride(4, 0));
      NS invalid{zero, ctx};
    } catch (const std::invalid_argument&) { rejected = true; }
    CHECK(rejected);

    rejected = false;
    try {
      // Both logical indices address the same element.
      Kokkos::LayoutStride overlap(2, 1, 2, 1);
      Kokkos::View<sunrealtype**, Kokkos::LayoutStride,
                   typename ExecSpace::memory_space> bad("bad", overlap);
      using NB = sundials::kokkos::VectorFor<decltype(bad)>;
      NB invalid{bad, ctx};
    } catch (const std::invalid_argument&) { rejected = true; }
    CHECK(rejected);

    // Managed views retain their allocation after the source handle goes away;
    // an unmanaged view is only borrowed and is never used after its owner.
    N1 retained{Kokkos::View<sunrealtype*, typename ExecSpace::memory_space>(
                  "retained", 4),
                ctx};
    auto retained_ptr = retained.View().data();
    CHECK(retained_ptr != nullptr);

    V2 shape_a("shape_a", 2, 4), shape_b("shape_b", 4, 2), shape_out("shape_out", 2, 4);
    N2 na{shape_a, ctx}, nb{shape_b, ctx}, no{shape_out, ctx};
    (void)SUNContext_GetLastError(ctx);
    N_VLinearSum(SUN_RCONST(1.0), na, SUN_RCONST(1.0), nb, no);
    CHECK(SUNContext_GetLastError(ctx) != SUN_SUCCESS);
    Kokkos::fence();
  }
  Kokkos::finalize();
  return 0;
}
