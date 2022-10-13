/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel McGreer, Cody Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file an NVector implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#ifndef _NVECTOR_KOKKOS_HPP
#define _NVECTOR_KOKKOS_HPP

#include <Kokkos_Core.hpp>
#include <memory>
#include <sundials/sundials_nvector.hpp>

namespace sundials {
namespace kokkos {

// Forward declaration
template<class ExecutionSpace, class MemorySpace>
class Vector;

// Get the Kokkos vector wrapped by an N_Vector
template<class ExecutionSpace,
         class MemorySpace = typename ExecutionSpace::memory_space>
inline Vector<ExecutionSpace, MemorySpace>* GetVec(N_Vector v)
{
  return static_cast<Vector<ExecutionSpace, MemorySpace>*>(v->content);
}

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

/* N_Vector API ops to implement.override */

SUNDIALS_INLINE
N_Vector_ID N_VGetVectorID_Kokkos(N_Vector v)
{
  return SUNDIALS_NVEC_KOKKOS;
}

template<class VectorType>
sunindextype N_VGetLength_Kokkos(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  return vec->Length();
}

template<class VectorType>
sunrealtype* N_VGetArrayPointer_Kokkos(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  return vec->HostView().data();
}

template<class VectorType>
sunrealtype* N_VGetDeviceArrayPointer_Kokkos(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  return vec->View().data();
}

template<class VectorType>
N_Vector N_VClone_Kokkos(N_Vector w)
{
  auto vec{static_cast<VectorType*>(w->content)};
  auto new_vec{new VectorType(*vec)}; // NOLINT
  return new_vec->Convert();
}

template<class VectorType>
void N_VDestroy_Kokkos(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  delete vec;
  return;
}

template<class VectorType>
void N_VPrint_Kokkos(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
}

template<class VectorType>
void N_VPrintFile_Kokkos(N_Vector v, FILE* outfile)
{
  auto vec{static_cast<VectorType*>(v->content)};
}

/* OPTIONAL local reduction kernels (no parallel communication) */

template<class VectorType>
sunrealtype N_VWSqrSumLocal_Kokkos(N_Vector x, N_Vector w)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto wvec{static_cast<VectorType*>(w->content)};
  auto wdata{wvec->View()};

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
      "N_VWSqrSumLocal", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) { update += (xdata(i) * wdata(i) * xdata(i) * wdata(i)); },
      gpu_result);

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VWSqrSumMaskLocal_Kokkos(N_Vector x, N_Vector w, N_Vector id)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto wvec{static_cast<VectorType*>(w->content)};
  auto wdata{wvec->View()};
  auto idvec{static_cast<VectorType*>(id->content)};
  auto iddata{idvec->View()};

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
      "N_VWSqrSumMaskLocal", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        if (iddata(i) > sunrealtype{0.0}) update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
      },
      gpu_result);

  return gpu_result;
}

/* standard vector operations */

template<class VectorType>
void N_VAbs_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};
  Kokkos::parallel_for(
      "N_VAbs", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = abs(xdata(i)); });
}

template<class VectorType>
void N_VAddConst_Kokkos(N_Vector x, sunrealtype b, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};
  Kokkos::parallel_for(
      "N_VAddConst", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = xdata(i) + b; });
}

template<class VectorType>
void N_VCompare_Kokkos(sunrealtype c, N_Vector x, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};
  Kokkos::parallel_for(
      "N_VCompare", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = abs(xdata(i)) >= c ? sunrealtype{1.0} : sunrealtype{0.0}; });
}

template<class VectorType>
void N_VConst_Kokkos(sunrealtype c, N_Vector z)
{
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};
  Kokkos::parallel_for(
      "N_VConst", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) { zdata(i) = c; });
}

template<class VectorType>
booleantype N_VConstrMask_Kokkos(N_Vector c, N_Vector x, N_Vector m)
{
  auto cvec{static_cast<VectorType*>(c->content)};
  auto cdata{cvec->View()};
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto mvec{static_cast<VectorType*>(m->content)};
  auto mdata{mvec->View()};

  sunrealtype sum{0.0};
  Kokkos::parallel_reduce(
      "N_VConstrMask", typename VectorType::range_policy(0, mvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        bool test = (abs(cdata(i)) > sunrealtype{1.5} && cdata(i) * xdata(i) <= sunrealtype{0.0}) ||
                    (abs(cdata(i)) > sunrealtype{0.5} && cdata(i) * xdata(i) < sunrealtype{0.0});
        mdata(i) = test ? sunrealtype{1.0} : sunrealtype{0.0};
        update += mdata(i);
      },
      sum);

  return (sum < sunrealtype{0.5});
}

template<class VectorType>
void N_VDiv_Kokkos(N_Vector x, N_Vector y, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto yvec{static_cast<VectorType*>(y->content)};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};
  Kokkos::parallel_for(
      "N_VDiv", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = xdata(i) / ydata(i); });
}

template<class VectorType>
sunrealtype N_VDotProd_Kokkos(N_Vector x, N_Vector y)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto yvec{static_cast<VectorType*>(y->content)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
      "N_VDotProd", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) { update += xdata(i) * ydata(i); }, gpu_result);

  return gpu_result;
}

template<class VectorType>
void N_VInv_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};

  Kokkos::parallel_for(
      "N_VInv", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = sunrealtype{1.0} / xdata(i); });
}

template<class VectorType>
booleantype N_VInvTest_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};

  sunrealtype minimum{0.0};
  Kokkos::parallel_reduce(
      "N_VInvTest", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        if (xdata(i) == sunrealtype{0.0}) {
          update += sunrealtype{1.0};
        }
        else {
          zdata(i) = sunrealtype{1.0} / xdata(i);
        }
      },
      minimum);

  return (minimum < sunrealtype{0.5});
}

template<class VectorType>
sunrealtype N_VL1Norm_Kokkos(N_Vector x)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
      "N_VL1Norm", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) { update += (abs(xdata(i))); }, gpu_result);

  return gpu_result;
}

template<class VectorType>
void N_VLinearSum_Kokkos(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto yvec{static_cast<VectorType*>(y->content)};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};

  Kokkos::parallel_for(
      "N_VLinearSum", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = a * xdata(i) + b * ydata(i); });
}

template<class VectorType>
sunrealtype N_VMaxNorm_Kokkos(N_Vector x)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
      "N_VMaxNorm", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        if (abs(xdata(i)) > update) update = abs(xdata(i));
      },
      Kokkos::Max<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VMin_Kokkos(N_Vector x)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};

  sunrealtype gpu_result{std::numeric_limits<sunrealtype>::max()};
  Kokkos::parallel_reduce(
      "N_VMin", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        if (xdata(i) < update) update = xdata(i);
      },
      Kokkos::Min<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VMinQuotient_Kokkos(N_Vector num, N_Vector denom)
{
  auto nvec{static_cast<VectorType*>(num->content)};
  auto ndata{nvec->View()};
  auto dvec{static_cast<VectorType*>(denom->content)};
  auto ddata{dvec->View()};

  sunrealtype gpu_result{std::numeric_limits<sunrealtype>::max()};

  Kokkos::parallel_reduce(
      "N_VMinQuotient", typename VectorType::range_policy(0, nvec->Length()),
      KOKKOS_LAMBDA(sunindextype i, sunrealtype & update) {
        if (ddata(i) != sunrealtype{0.0}) {
          if ((ndata(i) / ddata(i)) < update) update = ndata(i) / ddata(i);
        }
      },
      Kokkos::Min<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
void N_VProd_Kokkos(N_Vector x, N_Vector y, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto yvec{static_cast<VectorType*>(y->content)};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};

  Kokkos::parallel_for(
      "N_VProd", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = xdata(i) * ydata(i); });
}

template<class VectorType>
void N_VScale_Kokkos(sunrealtype c, N_Vector x, N_Vector z)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  auto xdata{xvec->View()};
  auto zvec{static_cast<VectorType*>(z->content)};
  auto zdata{zvec->View()};

  Kokkos::parallel_for(
      "N_VScale", typename VectorType::range_policy(0, xvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) { zdata(i) = c * xdata(i); });
}

template<class VectorType>
sunrealtype N_VWL2Norm_Kokkos(N_Vector x, N_Vector w)
{
  return std::sqrt(impl::N_VWSqrSumLocal_Kokkos<VectorType>(x, w));
}

template<class VectorType>
sunrealtype N_VWrmsNorm_Kokkos(N_Vector x, N_Vector w)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  return std::sqrt(impl::N_VWSqrSumLocal_Kokkos<VectorType>(x, w) / xvec->Length());
}

template<class VectorType>
sunrealtype N_VWrmsNormMask_Kokkos(N_Vector x, N_Vector w, N_Vector id)
{
  auto xvec{static_cast<VectorType*>(x->content)};
  return std::sqrt(impl::N_VWSqrSumMaskLocal_Kokkos<VectorType>(x, w, id) / xvec->Length());
}

/* fused vector operations */
template<class VectorType>
int N_VLinearCombination_Kokkos(int nvec, sunrealtype* c, N_Vector* X, N_Vector z)
{
  auto d_X{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_X[j] = static_cast<VectorType*>(X[j]->content)->View().data();

  auto zvec{static_cast<VectorType*>(z->content)};
  auto h_cview{Kokkos::View<sunrealtype*, Kokkos::HostSpace>(c, nvec)};
  auto cview{Kokkos::create_mirror_view_and_copy((typename VectorType::memory_space){}, h_cview)};
  auto zview{zvec->View()};
  Kokkos::parallel_for(
      "N_VLinearCombination", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        zview(i) = cview(0) * d_X[0][i];
        for (int j = 1; j < nvec; j++) {
          zview(i) += cview(j) * d_X[j][i];
        }
      });

  Kokkos::kokkos_free(d_X);

  return 0;
}

template<class VectorType>
int N_VScaleAddMulti_Kokkos(int nvec, sunrealtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  auto d_Y{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Y[j] = static_cast<VectorType*>(Y[j]->content)->View().data();

  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Z[j] = static_cast<VectorType*>(Z[j]->content)->View().data();

  auto xvec{static_cast<VectorType*>(x->content)};
  auto h_cview{Kokkos::View<sunrealtype*, Kokkos::HostSpace>(c, nvec)};
  auto cview{Kokkos::create_mirror_view_and_copy((typename VectorType::memory_space){}, h_cview)};
  auto xview{xvec->View()};
  Kokkos::parallel_for(
      "N_VScaleAddMulti", typename VectorType::range_policy(0, xvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++) d_Z[j][i] = cview(j) * xview(i) + d_Y[j][i];
      });

  Kokkos::kokkos_free(d_Y);
  Kokkos::kokkos_free(d_Z);

  return 0;
}

/* vector array operations */
template<class VectorType>
int N_VConstVectorArray_Kokkos(int nvec, sunrealtype c, N_Vector* Z)
{
  auto zvec{static_cast<VectorType*>(Z[0]->content)};
  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Z[j] = static_cast<VectorType*>(Z[j]->content)->View().data();

  Kokkos::parallel_for(
      "N_VConstVectorArray", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++) d_Z[j][i] = c;
      });

  Kokkos::kokkos_free(d_Z);

  return 0;
}

template<class VectorType>
int N_VLinearCombinationVectorArray_Kokkos(int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z)
{
  auto zvec{static_cast<VectorType*>(Z[0]->content)};

  auto d_X{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * nsum * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) {
    for (int k = 0; k < nsum; k++) {
      d_X[j * nsum + k] = static_cast<VectorType*>(X[k][j]->content)->View().data();
    }
  }

  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Z[j] = static_cast<VectorType*>(Z[j]->content)->View().data();

  auto h_cview{Kokkos::View<sunrealtype*, Kokkos::HostSpace>(c, nvec)};
  auto cview{Kokkos::create_mirror_view_and_copy((typename VectorType::memory_space){}, h_cview)};
  Kokkos::parallel_for(
      "N_VLinearCombinationVectorArray", typename VectorType::range_policy(0, zvec->Length()),
      KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++) {
          d_Z[j][i] = cview(0) * d_X[j * nsum][i];
          for (int k = 1; k < nsum; k++) {
            d_Z[j][i] += cview(k) * d_X[j * nsum + k][i];
          }
        }
      });

  Kokkos::kokkos_free(d_X);
  Kokkos::kokkos_free(d_Z);

  return 0;
}

template<class VectorType>
int N_VLinearSumVectorArray_Kokkos(int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  auto zvec{static_cast<VectorType*>(Z[0]->content)};

  auto d_X{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_X[j] = static_cast<VectorType*>(X[j]->content)->View().data();

  auto d_Y{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Y[j] = static_cast<VectorType*>(Y[j]->content)->View().data();

  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Z[j] = static_cast<VectorType*>(Z[j]->content)->View().data();

  Kokkos::parallel_for(
      "N_VLinearSumVectorArray", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++) d_Z[j][i] = a * d_X[j][i] + b * d_Y[j][i];
      });

  Kokkos::kokkos_free(d_X);
  Kokkos::kokkos_free(d_Y);
  Kokkos::kokkos_free(d_Z);

  return 0;
}

template<class VectorType>
int N_VScaleAddMultiVectorArray_Kokkos(int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  auto zvec{static_cast<VectorType*>(Z[0][0]->content)};

  auto d_X{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_X[j] = static_cast<VectorType*>(X[j]->content)->View().data();

  auto d_Y{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * nsum * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) {
    for (int k = 0; k < nsum; k++) {
      d_Y[j * nsum + k] = static_cast<VectorType*>(Y[k][j]->content)->View().data();
    }
  }

  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * nsum * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) {
    for (int k = 0; k < nsum; k++) {
      d_Z[j * nsum + k] = static_cast<VectorType*>(Z[k][j]->content)->View().data();
    }
  }

  auto h_aview{Kokkos::View<sunrealtype*, Kokkos::HostSpace>(a, nvec)};
  auto aview{Kokkos::create_mirror_view_and_copy((typename VectorType::memory_space){}, h_aview)};

  Kokkos::parallel_for(
      "N_VScaleAddMultiVectorArray", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++)
          for (int k = 0; k < nsum; k++) d_Z[j * nsum + k][i] = aview(k) * d_X[j][i] + d_Y[j * nsum + k][i];
      });

  Kokkos::kokkos_free(d_X);
  Kokkos::kokkos_free(d_Y);
  Kokkos::kokkos_free(d_Z);

  return 0;
}

template<class VectorType>
int N_VScaleVectorArray_Kokkos(int nvec, sunrealtype* c, N_Vector* X, N_Vector* Z)
{
  auto zvec{static_cast<VectorType*>(Z[0]->content)};

  auto d_X{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_X[j] = static_cast<VectorType*>(X[j]->content)->View().data();

  auto d_Z{static_cast<sunrealtype**>(Kokkos::kokkos_malloc(nvec * sizeof(sunrealtype*)))};
  for (int j = 0; j < nvec; j++) d_Z[j] = static_cast<VectorType*>(Z[j]->content)->View().data();

  auto h_cview{Kokkos::View<sunrealtype*, Kokkos::HostSpace>(c, nvec)};
  auto cview{Kokkos::create_mirror_view_and_copy((typename VectorType::memory_space){}, h_cview)};

  Kokkos::parallel_for(
      "N_VScaleVectorArray", typename VectorType::range_policy(0, zvec->Length()), KOKKOS_LAMBDA(sunindextype i) {
        for (int j = 0; j < nvec; j++) d_Z[j][i] = cview(j) * d_X[j][i];
      });

  return 0;
}

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombination_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Kokkos(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Kokkos(N_Vector v, booleantype tf);

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
         class MemorySpace = class ExecutionSpace::memory_space>
class Vector : public sundials::impl::BaseNvector,
               public sundials::ConvertibleTo<N_Vector>
{
public:
  using view_type      = Kokkos::View<sunrealtype*, MemorySpace>;
  using host_view_type = typename view_type::HostMirror;
  using memory_space   = MemorySpace;
  using exec_space     = typename MemorySpace::execution_space;
  using range_policy   = Kokkos::RangePolicy<exec_space>;

  Vector(sunindextype length, SUNContext sunctx)
      : length_(length), view_("device data", length),
        host_view_(Kokkos::create_mirror_view(view_)), sundials::impl::BaseNvector(sunctx)
  {
    initNvector();
  }

  Vector(sunindextype length, view_type view, SUNContext sunctx)
      : length_(length), view_(view), host_view_(Kokkos::create_mirror_view(view_)), sundials::impl::BaseNvector(sunctx)
  {
    initNvector();
  }

  Vector(sunindextype length, view_type view, host_view_type host_view, SUNContext sunctx)
      : length_(length), view_(view), host_view_(host_view), sundials::impl::BaseNvector(sunctx)
  {
    initNvector();
  }

  // Move constructor
  Vector(Vector&& that_vector) noexcept
      : length_(that_vector.length_), view_(std::move(that_vector.view_)),
        host_view_(std::move(that_vector.host_view_)), sundials::impl::BaseNvector(std::move(that_vector))
  {
    initNvector();
  }

  // Copy constructor
  Vector(const Vector& that_vector)
      : length_(that_vector.length_), view_("device data", length_),
        host_view_(), sundials::impl::BaseNvector(that_vector)
  {
    Kokkos::deep_copy(view_, that_vector.view_);
    host_view_ = Kokkos::create_mirror_view(view_);
    initNvector();
  }

  // Move assignment
  Vector& operator=(Vector&& rhs) noexcept
  {
    length_    = rhs.length_;
    view_      = std::move(rhs.view_);
    host_view_ = std::move(rhs.host_view_);

    sundials::impl::BaseNvector::operator=(std::move(rhs));

    return *this;
  }

  // Copy assignment
  Vector& operator=(const Vector& rhs)
  {
    length_ = rhs.length_;
    Kokkos::deep_copy(view_, rhs.view_);
    host_view_ = Kokkos::create_mirror_view(view_);

    sundials::impl::BaseNvector::operator=(rhs);

    return *this;
  }

  // Accessors

  sunindextype Length() const
  {
    return length_;
  }

  view_type View()
  {
    return view_;
  }

  host_view_type HostView()
  {
    return host_view_;
  }

  void EnableAllFusedOps(bool truefalse) {
    using this_type = Vector<ExecutionSpace, MemorySpace>;

    if (truefalse) {
      this->object_->ops->nvlinearcombination            = impl::N_VLinearCombination_Kokkos<this_type>;
      this->object_->ops->nvscaleaddmulti                = impl::N_VScaleAddMulti_Kokkos<this_type>;
      this->object_->ops->nvdotprodmulti                 = nullptr;
      this->object_->ops->nvlinearsumvectorarray         = impl::N_VLinearSumVectorArray_Kokkos<this_type>;
      this->object_->ops->nvscalevectorarray             = impl::N_VScaleVectorArray_Kokkos<this_type>;
      this->object_->ops->nvconstvectorarray             = impl::N_VConstVectorArray_Kokkos<this_type>;
      this->object_->ops->nvwrmsnormvectorarray          = nullptr;
      this->object_->ops->nvwrmsnormmaskvectorarray      = nullptr;
      this->object_->ops->nvscaleaddmultivectorarray     = impl::N_VScaleAddMultiVectorArray_Kokkos<this_type>;
      this->object_->ops->nvlinearcombinationvectorarray = impl::N_VLinearCombinationVectorArray_Kokkos<this_type>;
    } else {
      this->object_->ops->nvlinearcombination            = nullptr;
      this->object_->ops->nvscaleaddmulti                = nullptr;
      this->object_->ops->nvdotprodmulti                 = nullptr;
      this->object_->ops->nvlinearsumvectorarray         = nullptr;
      this->object_->ops->nvscalevectorarray             = nullptr;
      this->object_->ops->nvconstvectorarray             = nullptr;
      this->object_->ops->nvwrmsnormvectorarray          = nullptr;
      this->object_->ops->nvwrmsnormmaskvectorarray      = nullptr;
      this->object_->ops->nvscaleaddmultivectorarray     = nullptr;
      this->object_->ops->nvlinearcombinationvectorarray = nullptr;
    }
  }

  // Override ConvertibleTo operations

  operator N_Vector() override
  {
    return object_.get();
  }

  operator N_Vector() const override
  {
    return object_.get();
  }

  N_Vector Convert() override
  {
    return object_.get();
  }

  N_Vector Convert() const override
  {
    return object_.get();
  }

private:
  sunindextype length_;
  view_type view_;
  host_view_type host_view_;

  void initNvector()
  {
    using this_type = Vector<ExecutionSpace, MemorySpace>;

    this->object_->content = this;

    /* constructors, destructors, and utility operations */
    this->object_->ops->nvclone                 = impl::N_VClone_Kokkos<this_type>;
    this->object_->ops->nvdestroy               = impl::N_VDestroy_Kokkos<this_type>;
    this->object_->ops->nvgetarraypointer       = impl::N_VGetArrayPointer_Kokkos<this_type>;
    this->object_->ops->nvgetdevicearraypointer = impl::N_VGetDeviceArrayPointer_Kokkos<this_type>;
    this->object_->ops->nvgetlength             = impl::N_VGetLength_Kokkos<this_type>;
    this->object_->ops->nvgetvectorid           = impl::N_VGetVectorID_Kokkos;

    /* standard vector operations */
    this->object_->ops->nvabs          = impl::N_VAbs_Kokkos<this_type>;
    this->object_->ops->nvaddconst     = impl::N_VAddConst_Kokkos<this_type>;
    this->object_->ops->nvcompare      = impl::N_VCompare_Kokkos<this_type>;
    this->object_->ops->nvconst        = impl::N_VConst_Kokkos<this_type>;
    this->object_->ops->nvconstrmask   = impl::N_VConstrMask_Kokkos<this_type>;
    this->object_->ops->nvdiv          = impl::N_VDiv_Kokkos<this_type>;
    this->object_->ops->nvdotprod      = impl::N_VDotProd_Kokkos<this_type>;
    this->object_->ops->nvinv          = impl::N_VInv_Kokkos<this_type>;
    this->object_->ops->nvinvtest      = impl::N_VInvTest_Kokkos<this_type>;
    this->object_->ops->nvl1norm       = impl::N_VL1Norm_Kokkos<this_type>;
    this->object_->ops->nvlinearsum    = impl::N_VLinearSum_Kokkos<this_type>;
    this->object_->ops->nvmaxnorm      = impl::N_VMaxNorm_Kokkos<this_type>;
    this->object_->ops->nvmin          = impl::N_VMin_Kokkos<this_type>;
    this->object_->ops->nvminquotient  = impl::N_VMinQuotient_Kokkos<this_type>;
    this->object_->ops->nvprod         = impl::N_VProd_Kokkos<this_type>;
    this->object_->ops->nvscale        = impl::N_VScale_Kokkos<this_type>;
    this->object_->ops->nvwl2norm      = impl::N_VWL2Norm_Kokkos<this_type>;
    this->object_->ops->nvwrmsnorm     = impl::N_VWrmsNorm_Kokkos<this_type>;
    this->object_->ops->nvwrmsnormmask = impl::N_VWrmsNormMask_Kokkos<this_type>;

    /* local reduction operations */
    this->object_->ops->nvconstrmasklocal  = impl::N_VConstrMask_Kokkos<this_type>;
    this->object_->ops->nvdotprodlocal     = impl::N_VDotProd_Kokkos<this_type>;
    this->object_->ops->nvinvtestlocal     = impl::N_VInvTest_Kokkos<this_type>;
    this->object_->ops->nvl1normlocal      = impl::N_VL1Norm_Kokkos<this_type>;
    this->object_->ops->nvmaxnormlocal     = impl::N_VMaxNorm_Kokkos<this_type>;
    this->object_->ops->nvminlocal         = impl::N_VMin_Kokkos<this_type>;
    this->object_->ops->nvminquotientlocal = impl::N_VMinQuotient_Kokkos<this_type>;
    this->object_->ops->nvwsqrsumlocal     = impl::N_VWSqrSumLocal_Kokkos<this_type>;
    this->object_->ops->nvwsqrsummasklocal = impl::N_VWSqrSumMaskLocal_Kokkos<this_type>;
  }
};

/* Implementation specific */
template<class VectorType>
void CopyToDevice(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  CopyToDevice(*vec);
}

template<class VectorType>
void CopyFromDevice(N_Vector v)
{
  auto vec{static_cast<VectorType*>(v->content)};
  CopyFromDevice(*vec);
}

/* Implementation specific */
template<class VectorType>
void CopyToDevice(VectorType& v)
{
  Kokkos::deep_copy(v.View(), v.HostView());
}

template<class VectorType>
void CopyFromDevice(VectorType& v)
{
  Kokkos::deep_copy(v.HostView(), v.View());
}

} // namespace kokkos
} // namespace sundials

#endif
