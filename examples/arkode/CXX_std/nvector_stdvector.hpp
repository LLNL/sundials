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
 * module): a std::vector-based vector using the virtual sundials::NVectorBase.
 *
 * Usage:
 *
 *   // 0. Empty state (no N_Vector handle), e.g., for containers
 *   sundials::example::StdVector v;
 *
 *   // 1. Allocating
 *   sundials::example::StdVector v(100, sunctx);
 *
 *   // 2. Move from existing data
 *   std::vector<sunrealtype> my_data(100, 1.0);
 *   sundials::example::StdVector v(std::move(my_data), sunctx); // my_data is now empty
 *
 *   // 3. View existing data (no copy, no ownership transfer)
 *   std::vector<sunrealtype> my_data(100, 1.0);
 *   sundials::example::StdVector v(my_data, sunctx); // v reads/writes my_data directly
 *   // my_data must outlive v!
 *
 *   // 4. Clone (allocates a new owning vector; SUNDIALS frees it with
 *   //    N_VDestroy)
 *   sundials::example::StdVector* w = v.Clone(sunctx);
 *
 *   // 5. Access in a SUNDIALS callback
 *   int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
 *   {
 *     const auto* yvec = sundials::example::StdVector::FromNVectorConst(y);
 *     auto* ydotvec    = sundials::example::StdVector::FromNVector(ydot);
 *     ...
 *   }
 *
 * Copying a view copies the view (both alias the same user buffer). Assigning
 * an owning vector into a view makes the target owning.
 *
 * Do not move out of an object that SUNDIALS still holds (e.g., one obtained
 * with FromNVector from a clone, or one whose N_Vector was given to an
 * integrator): the handle would then refer to a moved-from object. A moved-from
 * owning StdVector reports length 0 and may be destroyed or assigned to.
 *
 * Operations that are overridden below are installed in the N_Vector; all
 * others are left NULL (see sundials_nvector_base.hpp).
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_STDVECTOR_HPP
#define SUNDIALS_NVECTOR_STDVECTOR_HPP

#include <sundials/sundials_nvector_base.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace sundials {
namespace example {

class StdVector final : public NVectorBase<StdVector>
{
public:
  // ---------------------------------------------------------------------------
  // Construction
  // ---------------------------------------------------------------------------

  // Empty state: no N_Vector handle and no data.
  StdVector() = default;

  // Allocating constructor. Creates and owns a vector.
  StdVector(sunindextype num_elements, SUNContext sunctx)
    : NVectorBase(sunctx), owned_(static_cast<std::size_t>(num_elements))
  {}

  // Move constructor. Takes ownership of an existing std::vector. The source
  // is left in a moved-from state (typically empty).
  StdVector(std::vector<sunrealtype>&& values, SUNContext sunctx)
    : NVectorBase(sunctx),
      owned_(std::move(values))
  {}

  // View constructor. Wraps an existing std::vector WITHOUT copying or taking
  // ownership. The caller must ensure the source outlives this StdVector.
  // Reads and writes go directly to the source vector's storage.
  StdVector(std::vector<sunrealtype>& values, SUNContext sunctx)
    : NVectorBase(sunctx),
      view_(values.data()),
      view_len_(static_cast<sunindextype>(values.size()))
  {}

  // Clones always allocate fresh owning storage. SUNDIALS takes ownership of
  // the clone and frees it via N_VDestroy.
  StdVector* Clone(SUNContext sunctx) const override
  {
    return new StdVector(length(), sunctx);
  }

  // ---------------------------------------------------------------------------
  // Utilities
  // ---------------------------------------------------------------------------

  sunindextype GetLength() const override { return length(); }

  // Returns a pointer to the underlying data. Valid for both owning and view
  // modes.
  sunrealtype* data() { return ptr(); }

  const sunrealtype* data() const { return ptr(); }

  // Convenience operator for element access.
  sunrealtype& operator[](sunindextype i) { return ptr()[i]; }

  const sunrealtype& operator[](sunindextype i) const { return ptr()[i]; }

  // Returns true if this StdVector owns its data (allocated or moved-in).
  // Returns false for views.
  bool owning() const { return view_ == nullptr; }

  sunrealtype* GetArrayPointer() override { return ptr(); }

  // ---------------------------------------------------------------------------
  // Standard vector operations
  // ---------------------------------------------------------------------------

  void LinearSum(sunrealtype a, const StdVector& x, sunrealtype b,
                 const StdVector& y) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    const sunrealtype* ydata = y.ptr();
    for (sunindextype i = 0; i < length(); ++i)
      zdata[i] = a * xdata[i] + b * ydata[i];
  }

  void Const(sunrealtype c) override { std::fill(ptr(), ptr() + length(), c); }

  void Prod(const StdVector& x, const StdVector& y) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    const sunrealtype* ydata = y.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = xdata[i] * ydata[i];
  }

  void Div(const StdVector& x, const StdVector& y) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    const sunrealtype* ydata = y.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = xdata[i] / ydata[i];
  }

  void Scale(sunrealtype c, const StdVector& x) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = c * xdata[i];
  }

  void Abs(const StdVector& x) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = std::abs(xdata[i]);
  }

  void Inv(const StdVector& x) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = SUN_RCONST(1.0) / xdata[i];
  }

  void AddConst(const StdVector& x, sunrealtype b) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i) zdata[i] = xdata[i] + b;
  }

  sunrealtype DotProd(const StdVector& y) const override
  {
    const sunrealtype* xdata = ptr();
    const sunrealtype* ydata = y.ptr();
    sunrealtype sum          = 0;
    for (sunindextype i = 0; i < length(); ++i) sum += xdata[i] * ydata[i];
    return sum;
  }

  sunrealtype MaxNorm() const override
  {
    const sunrealtype* xdata = ptr();
    sunrealtype max_abs      = 0;
    for (sunindextype i = 0; i < length(); ++i)
      max_abs = std::max(max_abs, std::abs(xdata[i]));
    return max_abs;
  }

  sunrealtype WrmsNorm(const StdVector& w) const override
  {
    return std::sqrt(WSqrSumLocal(w) / static_cast<sunrealtype>(length()));
  }

  sunrealtype WrmsNormMask(const StdVector& w, const StdVector& id) const override
  {
    return std::sqrt(WSqrSumMaskLocal(w, id) / static_cast<sunrealtype>(length()));
  }

  sunrealtype Min() const override
  {
    const sunrealtype* xdata = ptr();
    sunrealtype min_val      = std::numeric_limits<sunrealtype>::max();
    for (sunindextype i = 0; i < length(); ++i) min_val = std::min(min_val, xdata[i]);
    return min_val;
  }

  sunrealtype WL2Norm(const StdVector& w) const override
  {
    return std::sqrt(WSqrSumLocal(w));
  }

  sunrealtype L1Norm() const override
  {
    const sunrealtype* xdata = ptr();
    sunrealtype sum          = 0;
    for (sunindextype i = 0; i < length(); ++i) sum += std::abs(xdata[i]);
    return sum;
  }

  void Compare(sunrealtype c, const StdVector& x) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i)
      zdata[i] = std::abs(xdata[i]) >= c ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
  }

  sunbooleantype InvTest(const StdVector& x) override
  {
    sunrealtype* zdata       = ptr();
    const sunrealtype* xdata = x.ptr();
    for (sunindextype i = 0; i < length(); ++i)
    {
      if (xdata[i] == SUN_RCONST(0.0)) return SUNFALSE;
      zdata[i] = SUN_RCONST(1.0) / xdata[i];
    }
    return SUNTRUE;
  }

  sunbooleantype ConstrMask(const StdVector& c, const StdVector& x) override
  {
    sunrealtype* mdata       = ptr();
    const sunrealtype* cdata = c.ptr();
    const sunrealtype* xdata = x.ptr();
    sunrealtype sum          = 0;
    for (sunindextype i = 0; i < length(); ++i)
    {
      bool test = (std::abs(cdata[i]) > SUN_RCONST(1.5) && cdata[i] * xdata[i] <= SUN_RCONST(0.0)) ||
                  (std::abs(cdata[i]) > SUN_RCONST(0.5) && cdata[i] * xdata[i] < SUN_RCONST(0.0));
      mdata[i]  = test ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
      sum += mdata[i];
    }
    return (sum < SUN_RCONST(0.5));
  }

  sunrealtype MinQuotient(const StdVector& denom) const override
  {
    const sunrealtype* numdata   = ptr();
    const sunrealtype* denomdata = denom.ptr();
    sunrealtype min_quot         = std::numeric_limits<sunrealtype>::max();
    for (sunindextype i = 0; i < length(); ++i)
      if (denomdata[i] != SUN_RCONST(0.0))
        min_quot = std::min(min_quot, numdata[i] / denomdata[i]);
    return min_quot;
  }

  // ---------------------------------------------------------------------------
  // Fused ops
  //
  // Outputs may alias inputs (see the ALIASING section of
  // sundials_nvector_base.hpp), so every loop reads all inputs at index i
  // before writing the output at index i.
  // ---------------------------------------------------------------------------

  SUNErrCode LinearCombination(ConstReals c, ConstVectors X) override
  {
    sunrealtype* zdata = ptr();
    for (sunindextype i = 0; i < length(); ++i)
    {
      sunrealtype sum = 0;
      for (std::size_t j = 0; j < c.size(); ++j) sum += c[j] * X[j].ptr()[i];
      zdata[i] = sum;
    }
    return SUN_SUCCESS;
  }

  SUNErrCode ScaleAddMulti(ConstReals a, ConstVectors Y, Vectors Z) const override
  {
    const sunrealtype* xdata = ptr();
    for (sunindextype i = 0; i < length(); ++i)
    {
      const sunrealtype x_i = xdata[i];
      for (std::size_t j = 0; j < Y.size(); ++j)
        Z[j].ptr()[i] = a[j] * x_i + Y[j].ptr()[i];
    }
    return SUN_SUCCESS;
  }

  SUNErrCode DotProdMulti(ConstVectors Y, Reals dotprods) const override
  {
    return DotProdMultiLocal(Y, dotprods);
  }

  // ---------------------------------------------------------------------------
  // Vector array ops
  // ---------------------------------------------------------------------------

  SUNErrCode LinearSumVectorArray(sunrealtype a, ConstVectors X, sunrealtype b,
                                  ConstVectors Y, Vectors Z) override
  {
    for (std::size_t j = 0; j < X.size(); ++j) Z[j].LinearSum(a, X[j], b, Y[j]);
    return SUN_SUCCESS;
  }

  SUNErrCode ScaleVectorArray(ConstReals c, ConstVectors X, Vectors Z) override
  {
    for (std::size_t j = 0; j < X.size(); ++j) Z[j].Scale(c[j], X[j]);
    return SUN_SUCCESS;
  }

  SUNErrCode ConstVectorArray(sunrealtype c, Vectors Z) override
  {
    for (std::size_t j = 0; j < Z.size(); ++j) Z[j].Const(c);
    return SUN_SUCCESS;
  }

  SUNErrCode WrmsNormVectorArray(ConstVectors X, ConstVectors W,
                                 Reals nrm) const override
  {
    for (std::size_t j = 0; j < X.size(); ++j) nrm[j] = X[j].WrmsNorm(W[j]);
    return SUN_SUCCESS;
  }

  SUNErrCode WrmsNormMaskVectorArray(ConstVectors X, ConstVectors W,
                                     const StdVector& id,
                                     Reals nrm) const override
  {
    for (std::size_t j = 0; j < X.size(); ++j)
      nrm[j] = X[j].WrmsNormMask(W[j], id);
    return SUN_SUCCESS;
  }

  SUNErrCode ScaleAddMultiVectorArray(ConstReals a, ConstVectors X,
                                      ConstVectors2D Y,
                                      Vectors2D Z) const override
  {
    for (std::size_t j = 0; j < Y.rows(); ++j)
      for (std::size_t i = 0; i < X.size(); ++i)
        Z[j][i].LinearSum(a[j], X[i], SUN_RCONST(1.0), Y[j][i]);
    return SUN_SUCCESS;
  }

  SUNErrCode LinearCombinationVectorArray(ConstReals c, ConstVectors2D X,
                                          Vectors Z) override
  {
    for (std::size_t i = 0; i < Z.size(); ++i)
    {
      sunrealtype* zdata = Z[i].ptr();
      for (sunindextype k = 0; k < Z[i].length(); ++k)
      {
        sunrealtype sum = 0;
        for (std::size_t j = 0; j < c.size(); ++j) sum += c[j] * X[j][i].ptr()[k];
        zdata[k] = sum;
      }
    }
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Local reduction kernels
  // ---------------------------------------------------------------------------

  sunrealtype WSqrSumLocal(const StdVector& w) const override
  {
    const sunrealtype* xdata = ptr();
    const sunrealtype* wdata = w.ptr();
    sunrealtype sum          = 0;
    for (sunindextype i = 0; i < length(); ++i)
    {
      const sunrealtype prod = xdata[i] * wdata[i];
      sum += prod * prod;
    }
    return sum;
  }

  sunrealtype WSqrSumMaskLocal(const StdVector& w,
                               const StdVector& id) const override
  {
    const sunrealtype* xdata  = ptr();
    const sunrealtype* wdata  = w.ptr();
    const sunrealtype* iddata = id.ptr();
    sunrealtype sum           = 0;
    for (sunindextype i = 0; i < length(); ++i)
      if (iddata[i] > SUN_RCONST(0.0))
      {
        const sunrealtype prod = xdata[i] * wdata[i];
        sum += prod * prod;
      }
    return sum;
  }

  // ---------------------------------------------------------------------------
  // Single buffer reduction
  // ---------------------------------------------------------------------------

  SUNErrCode DotProdMultiLocal(ConstVectors Y, Reals dotprods) const override
  {
    const sunrealtype* xdata = ptr();
    for (std::size_t j = 0; j < Y.size(); ++j)
    {
      const sunrealtype* ydata = Y[j].ptr();
      sunrealtype sum          = 0;
      for (sunindextype i = 0; i < length(); ++i) sum += xdata[i] * ydata[i];
      dotprods[j] = sum;
    }
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Buffer Pack/Unpack
  // ---------------------------------------------------------------------------

  SUNErrCode BufSize(sunindextype* size) const override
  {
    *size = length() * static_cast<sunindextype>(sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode BufPack(void* buf) const override
  {
    std::memcpy(buf, ptr(), static_cast<std::size_t>(length()) * sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode BufUnpack(void* buf) override
  {
    std::memcpy(ptr(), buf, static_cast<std::size_t>(length()) * sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Print
  // ---------------------------------------------------------------------------

  void Print() const override { PrintFile(stdout); }

  void PrintFile(FILE* outfile) const override
  {
    const sunrealtype* xdata = ptr();
    for (sunindextype i = 0; i < length(); ++i)
      std::fprintf(outfile, "  data[%ld] = " SUN_FORMAT_G "\n", static_cast<long>(i), xdata[i]);
  }

private:
  // Number of elements. Derived from the storage (not stored) for owning
  // vectors, so a moved-from vector reports the length of its (empty) storage.
  sunindextype length() const
  {
    return view_ ? view_len_ : static_cast<sunindextype>(owned_.size());
  }

  // Active data pointer: the user's buffer in view mode, else the owned one.
  sunrealtype* ptr() { return view_ ? view_ : owned_.data(); }

  const sunrealtype* ptr() const { return view_ ? view_ : owned_.data(); }

  std::vector<sunrealtype> owned_; // Storage when owning (empty for views)
  sunrealtype* view_ = nullptr;    // Non-null only in view mode
  sunindextype view_len_ = 0;      // Number of elements in view mode
};

static_assert(std::is_nothrow_move_constructible_v<StdVector>,
              "StdVector must be nothrow move constructible so std::vector "
              "growth keeps N_Vector handles valid");

} // namespace example
} // namespace sundials

#endif // SUNDIALS_NVECTOR_STDVECTOR_HPP
