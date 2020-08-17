#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_kokkos.h>

#include "sundials_debug.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

// Static constants
static constexpr sunindextype zeroIdx = 0;

// Helpful macros
#define NVEC_KOKKOS_CONTENT(x) ((N_VectorContent_Kokkos)(x->content))
#define NVEC_KOKKOS_MEMSIZE(x) (NVEC_KOKKOS_CONTENT(x)->length * sizeof(realtype))

static void AllocateData(N_Vector v);

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Kokkos(N_Vector v)
{
  return SUNDIALS_NVEC_KOKKOS;
}

N_Vector N_VNewEmpty_Kokkos()
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Kokkos;
  v->ops->nvclone           = N_VClone_Kokkos;
  v->ops->nvcloneempty      = N_VCloneEmpty_Kokkos;
  v->ops->nvdestroy         = N_VDestroy_Kokkos;
  v->ops->nvspace           = N_VSpace_Kokkos;
  v->ops->nvgetlength       = N_VGetLength_Kokkos;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Kokkos;
  v->ops->nvconst        = N_VConst_Kokkos;
  v->ops->nvprod         = N_VProd_Kokkos;
  v->ops->nvdiv          = N_VDiv_Kokkos;
  v->ops->nvscale        = N_VScale_Kokkos;
  v->ops->nvabs          = N_VAbs_Kokkos;
  v->ops->nvinv          = N_VInv_Kokkos;
  v->ops->nvaddconst     = N_VAddConst_Kokkos;
  v->ops->nvdotprod      = N_VDotProd_Kokkos;
  v->ops->nvmaxnorm      = N_VMaxNorm_Kokkos;
  v->ops->nvmin          = N_VMin_Kokkos;
  v->ops->nvl1norm       = N_VL1Norm_Kokkos;
  v->ops->nvinvtest      = N_VInvTest_Kokkos;
  v->ops->nvconstrmask   = N_VConstrMask_Kokkos;
  v->ops->nvminquotient  = N_VMinQuotient_Kokkos;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Kokkos;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Kokkos;
  v->ops->nvwl2norm      = N_VWL2Norm_Kokkos;
  v->ops->nvcompare      = N_VCompare_Kokkos;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Kokkos;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Kokkos;
  v->ops->nvdotprodlocal     = N_VDotProd_Kokkos;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Kokkos;
  v->ops->nvminlocal         = N_VMin_Kokkos;
  v->ops->nvl1normlocal      = N_VL1Norm_Kokkos;
  v->ops->nvinvtestlocal     = N_VInvTest_Kokkos;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Kokkos;
  v->ops->nvminquotientlocal = N_VMinQuotient_Kokkos;

  v->content = (N_VectorContent_Kokkos) malloc(sizeof(_N_VectorContent_Kokkos));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }
}

N_Vector N_VNew_Kokkos(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Kokkos();
  if (v == NULL) return(NULL);

  NVEC_KOKKOS_CONTENT(v)->length = length;

  AllocateData(v);

  return(v);
}

N_Vector N_VMake_Kokkos(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Kokkos();
  if (v == NULL) return(NULL);

  NVEC_KOKKOS_CONTENT(v)->length      = length;
  // Create "unmanaged views" for data by passing data pointer into view contructor
  NVEC_KOKKOS_CONTENT(v)->host_data   = HostArrayView(h_vdata, length);
  NVEC_KOKKOS_CONTENT(v)->device_data = DeviceArrayView(d_vdata, length);

  return(v);
}

/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 */
sunindextype N_VGetLength_Kokkos(N_Vector v)
{
  return NVEC_KOKKOS_CONTENT(v)->length;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

// TODO, this currently returns actual data pointer
// Consider adding functions to return views?
// Same with get device array
realtype *N_VGetHostArrayPointer_Kokkos(N_Vector x)
{
  return NVEC_KOKKOS_CONTENT(x)->host_data.data();
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Kokkos(N_Vector x)
{
  return NVEC_KOKKOS_CONTENT(x)->device_data.data();;
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Kokkos(N_Vector x)
{
  Kokkos::deep_copy(NVEC_KOKKOS_CONTENT(x)->device_data,
                    NVEC_KOKKOS_CONTENT(x)->host_data);
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Kokkos(N_Vector x)
{

  Kokkos::deep_copy(NVEC_KOKKOS_CONTENT(x)->host_data,
                    NVEC_KOKKOS_CONTENT(x)->device_data);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Kokkos(N_Vector X)
{
  N_VPrintFile_Kokkos(X, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Kokkos(N_Vector X, FILE *outfile)
{
  const realtype *xd = NVEC_KOKKOS_CONTENT(X)->host_data.data();
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  sunindextype i;

  for (i = 0; i < N; ++i) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd[i]);
#else
    fprintf(outfile, "%11.8g\n", xd[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Kokkos(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Kokkos();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_KOKKOS_CONTENT(v)->length   = NVEC_KOKKOS_CONTENT(w)->length;

  return(v);
}

N_Vector N_VClone_Kokkos(N_Vector w)
{
  N_Vector v;
  v = NULL;
  v = N_VCloneEmpty_Kokkos(w);
  if (v == NULL) return(NULL);

  AllocateData(v);

  return(v);
}


void N_VDestroy_Kokkos(N_Vector v)
{
  if (v == NULL) return;

  N_VectorContent_Kokkos vc = NVEC_KOKKOS_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

   /*
   By assigning the views to default constructed views the refference
   count of the orignal view goes to zero and it automatically deallocates
   the data, unless the vector doesn't own the data in which it is a unmanaged
   view and the data will not be deallocated
   */
  vc->device_data = DeviceArrayView();
  vc->host_data = HostArrayView();

  /* free content struct */
  free(vc);
  v->content = NULL;

  /* free ops */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* free vector */
  free(v);
  v = NULL;

  return;
}


void N_VSpace_Kokkos(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_KOKKOS_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Kokkos(realtype c, N_Vector Z)
{
  const sunindextype N = NVEC_KOKKOS_CONTENT(Z)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VConst", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = c;
    }
  );

}

void N_VLinearSum_Kokkos(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;


  Kokkos::parallel_for("N_VLinearSum", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = a*xdata(i) + b*ydata(i);
    }
  );
}

void N_VProd_Kokkos(N_Vector X, N_Vector Y, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VProd", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = xdata(i) * ydata(i);
    }
  );
}


void N_VDiv_Kokkos(N_Vector X, N_Vector Y, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VDiv", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = xdata(i) / ydata(i);
    }
  );
}

void N_VScale_Kokkos(realtype c, N_Vector X, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VScale", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = c * xdata(i);
    }
  );
}

void N_VAbs_Kokkos(N_Vector X, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VAbs", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = abs(xdata(i));
    }
  );
}

void N_VInv_Kokkos(N_Vector X, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VInv", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = ONE / xdata(i);
    }
  );
}

void N_VAddConst_Kokkos(N_Vector X, realtype b, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VAddConst", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = xdata(i) + b;
    }
  );
}

realtype N_VDotProd_Kokkos(N_Vector X, N_Vector Y)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto ydata = NVEC_KOKKOS_CONTENT(Y)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VDotProd", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      update += xdata(i) * ydata(i);
    }, gpu_result);

  return (static_cast<realtype>(gpu_result));
}

realtype N_VMaxNorm_Kokkos(N_Vector X)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VMaxNorm", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if(abs(xdata(i)) > update) update = abs(xdata(i));
    }, Kokkos::Max<realtype>(gpu_result));

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWSqrSumLocal_Kokkos(N_Vector X, N_Vector W)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto wdata = NVEC_KOKKOS_CONTENT(W)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VWSqrSumLocal", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
    update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    }, gpu_result);

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNorm_Kokkos(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Kokkos(X, W);
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VWSqrSumMaskLocal_Kokkos(N_Vector X, N_Vector W, N_Vector ID)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const auto wdata = NVEC_KOKKOS_CONTENT(W)->device_data;
  const auto iddata = NVEC_KOKKOS_CONTENT(ID)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VWSqrSumMaskLocal", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (iddata(i) > ZERO)
        update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    }, gpu_result );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNormMask_Kokkos(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype sum = N_VWSqrSumMaskLocal_Kokkos(X, W, ID);
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VMin_Kokkos(N_Vector X)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = std::numeric_limits<realtype>::max();

  Kokkos::parallel_reduce( "N_VMin", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if(xdata(i) < update) update = xdata(i);
    }, Kokkos::Min<realtype>(gpu_result));

  return (static_cast<realtype>(gpu_result));

}

realtype N_VWL2Norm_Kokkos(N_Vector X, N_Vector W)
{
  return std::sqrt(N_VWSqrSumLocal_Kokkos(X, W));
}

realtype N_VL1Norm_Kokkos(N_Vector X)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VL1Norm", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      update += (abs(xdata(i)));
    }, gpu_result );

  return (static_cast<realtype>(gpu_result));
}

void N_VCompare_Kokkos(realtype c, N_Vector X, N_Vector Z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(X)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(X)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(Z)->device_data;

  Kokkos::parallel_for("N_VCompare", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      zdata(i) = abs(xdata(i)) >= c ? ONE : ZERO;
    }
  );

}

booleantype N_VInvTest_Kokkos(N_Vector x, N_Vector z)
{
  const auto xdata = NVEC_KOKKOS_CONTENT(x)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(x)->length;
  auto zdata = NVEC_KOKKOS_CONTENT(z)->device_data;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce( "N_VInvTest", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (xdata(i) == ZERO) {
        update += ONE;
      } else {
        zdata(i) = ONE/xdata(i);
      }
    }, gpu_result);

  realtype minimum = static_cast<realtype>(gpu_result);
  return (minimum < HALF);
}

booleantype N_VConstrMask_Kokkos(N_Vector c, N_Vector x, N_Vector m)
{
  const auto cdata = NVEC_KOKKOS_CONTENT(c)->device_data;
  const auto xdata = NVEC_KOKKOS_CONTENT(x)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(x)->length;
  auto mdata = NVEC_KOKKOS_CONTENT(m)->device_data;

  realtype gpu_result = 0.0;

  Kokkos::parallel_reduce("N_VConstrMask", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
    /*  bool test = (abs(cdata(i)) > ONEPT5 && cdata(i)*xdata(i) <= ZERO) ||
                  (abs(cdata(i)) > HALF   && cdata(i)*xdata(i) <  ZERO);
      mdata(i) = test ? ONE : ZERO;
      */
      update += mdata(i);
    }, gpu_result);

  realtype sum = static_cast<realtype>(gpu_result);
  return(sum < HALF);
}


realtype N_VMinQuotient_Kokkos(N_Vector num, N_Vector denom)
{
  const auto ndata = NVEC_KOKKOS_CONTENT(num)->device_data;
  const auto ddata = NVEC_KOKKOS_CONTENT(denom)->device_data;
  const sunindextype N = NVEC_KOKKOS_CONTENT(num)->length;

  realtype gpu_result = std::numeric_limits<realtype>::max();

  Kokkos::parallel_reduce( "N_VMinQuotient", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i, realtype &update) {
      if (ddata(i) != ZERO){
        if((ndata(i)/ddata(i)) < update) update = ndata(i)/ddata(i);
      }
    }, Kokkos::Min<realtype>(gpu_result));

  return (static_cast<realtype>(gpu_result));
}

/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

//TODO All functions take in realtype pointers
// maybe some implementation could take in Views?
int N_VLinearCombination_Kokkos(int nvec, realtype* c, N_Vector* X, N_Vector z)
{

  sunindextype N = NVEC_KOKKOS_CONTENT(z)->length;
  auto d_zd = NVEC_KOKKOS_CONTENT(z)->device_data;

  // Make c into a view, create a device view d_c, then deep copy c to it
  HostArrayView h_c(c, nvec);
  DeviceArrayView d_c("d_c", nvec);
  Kokkos::deep_copy(d_c, h_c);

  // Create View of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Xd("h_Xd", nvec);
  for (int j=0; j<nvec; j++)
    h_Xd(j) = NVEC_KOKKOS_CONTENT(X[j])->device_data;

  // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Xd("d_Xd", nvec);
  Kokkos::deep_copy(d_Xd, h_Xd);

  Kokkos::parallel_for("N_VLinearCombination", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      d_zd(i) = d_c(0) * d_Xd(0)(i);
      for (int j=1; j<nvec; j++){
        d_zd(i) += d_c(j) * d_Xd(j)(i);
      }
    }
  );

  return(0);
}

int N_VScaleAddMulti_Kokkos(int nvec, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(x)->length;
  auto d_xd = NVEC_KOKKOS_CONTENT(x)->device_data;

  // Make c into a view, create a device view d_c, then deep copy c to it,
  HostArrayView h_c(c, nvec);
  DeviceArrayView d_c("d_c", nvec);
  Kokkos::deep_copy(d_c, h_c);

  // Create array of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Yd("h_Yd", nvec);
  for (int j=0; j<nvec; j++)
    h_Yd(j) = NVEC_KOKKOS_CONTENT(Y[j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nvec);
  for (int j=0; j<nvec; j++)
    h_Zd(j) = NVEC_KOKKOS_CONTENT(Z[j])->device_data;

  // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Yd("d_Yd", nvec);
  Kokkos::deep_copy(d_Yd, h_Yd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for("N_VScaleAddMulti", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++)
        d_Zd(j)(i) = d_c(j) * d_xd(i) + d_Yd(j)(i);
    }
  );

  return(0);
}

/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Kokkos(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(Z[0])->length;

  // Create array of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Xd("h_Xd", nvec);
  for (int j=0; j<nvec; j++)
    h_Xd(j) = NVEC_KOKKOS_CONTENT(X[j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Yd("h_Yd", nvec);
  for (int j=0; j<nvec; j++)
    h_Yd(j) = NVEC_KOKKOS_CONTENT(Y[j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nvec);
  for (int j=0; j<nvec; j++)
    h_Zd(j) = NVEC_KOKKOS_CONTENT(Z[j])->device_data;

  // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Xd("d_Xd", nvec);
  Kokkos::deep_copy(d_Xd, h_Xd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Yd("d_Yd", nvec);
  Kokkos::deep_copy(d_Yd, h_Yd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for(" N_VLinearSumVectorArray", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++)
        d_Zd(j)(i) = a * d_Xd(j)(i) + b * d_Yd(j)(i);
    }
  );

  return(0);
}


int N_VScaleVectorArray_Kokkos(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(Z[0])->length;

  // Make c into a view, create a device view d_c, then deep copy c to it,
  HostArrayView h_c(c, nvec);
  DeviceArrayView d_c("d_c", nvec);
  Kokkos::deep_copy(d_c, h_c);

  // Create array of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Xd("h_Xd", nvec);
  for (int j=0; j<nvec; j++)
    h_Xd(j) = NVEC_KOKKOS_CONTENT(X[j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nvec);
  for (int j=0; j<nvec; j++)
    h_Zd(j) = NVEC_KOKKOS_CONTENT(Z[j])->device_data;

  // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Xd("d_Xd", nvec);
  Kokkos::deep_copy(d_Xd, h_Xd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for("N_VScaleVectorArray", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++)
        d_Zd(j)(i) = d_c(j) * d_Xd(j)(i);
    }
  );

  return(0);
}

int N_VConstVectorArray_Kokkos(int nvec, realtype c, N_Vector* Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(Z[0])->length;

  // Create array of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nvec);
  for (int j=0; j<nvec; j++)
    h_Zd(j) = NVEC_KOKKOS_CONTENT(Z[j])->device_data;

   // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for("N_VConstVectorArray", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++)
        d_Zd(j)(i) = c;
      }
    );

  return(0);
}

int N_VScaleAddMultiVectorArray_Kokkos(int nvec, int nsum, realtype* c,
                                       N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(X[0])->length;

  // Make c into a view, create a device view d_c, then deep copy c to it,
  HostArrayView h_c(c, nsum);
  DeviceArrayView d_c("d_c", nsum);
  Kokkos::deep_copy(d_c, h_c);

  // Create array of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Xd("h_Xd", nvec);
  for (int j=0; j<nvec; j++)
    h_Xd(j) = NVEC_KOKKOS_CONTENT(X[j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Yd("h_Yd", nsum*nvec);
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Yd(j*nsum+k) = NVEC_KOKKOS_CONTENT(Y[k][j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nsum*nvec);
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Zd(j*nsum+k) = NVEC_KOKKOS_CONTENT(Z[k][j])->device_data;

   // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Xd("d_Xd", nvec);
  Kokkos::deep_copy(d_Xd, h_Xd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Yd("d_Yd", nsum*nvec);
  Kokkos::deep_copy(d_Yd, h_Yd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nsum*nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for("N_VScaleAddMultiVectorArray", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++)
        for (int k=0; k<nsum; k++)
          d_Zd(j*nsum+k)(i) = d_c(k) * d_Xd(j)(i) + d_Yd(j*nsum+k)(i);
    }
  );

  return(0);
}


int N_VLinearCombinationVectorArray_Kokkos(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  sunindextype N = NVEC_KOKKOS_CONTENT(Z[0])->length;

  // Make c into a view, create a device view d_c, then deep copy c to it,
  HostArrayView h_c(c, nsum);
  DeviceArrayView d_c("d_c", nsum);
  Kokkos::deep_copy(d_c, h_c);

  // Create View of device views on host
  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Xd("h_Xd", nsum*nvec);
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Xd(j*nsum+k) = NVEC_KOKKOS_CONTENT(X[k][j])->device_data;

  Kokkos::View<DeviceArrayView*, Kokkos::HostSpace> h_Zd("h_Zd", nvec);
  for (int j=0; j<nvec; j++)
    h_Zd(j) = NVEC_KOKKOS_CONTENT(Z[j])->device_data;

  // Copy host view to a device view
  Kokkos::View<DeviceArrayView*, MemSpace> d_Xd("d_Xd", nsum*nvec);
  Kokkos::deep_copy(d_Xd, h_Xd);
  Kokkos::View<DeviceArrayView*, MemSpace> d_Zd("d_Zd", nvec);
  Kokkos::deep_copy(d_Zd, h_Zd);

  Kokkos::parallel_for("N_VLinearCombinationVectorArray", range_policy(zeroIdx, N),
    KOKKOS_LAMBDA(sunindextype i){
      for (int j=0; j<nvec; j++) {
        d_Zd(j)(i) = d_c(0) * d_Xd(j*nsum)(i);
        for (int k=1; k<nsum; k++) {
          d_Zd(j)(i) += d_c(k) * d_Xd(j*nsum+k)(i);
        }
      }
    }
  );

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Kokkos;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Kokkos;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Kokkos;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Kokkos;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Kokkos;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Kokkos;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Kokkos;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Kokkos;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Kokkos;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Kokkos;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Kokkos;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Kokkos;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Kokkos;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Kokkos(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Kokkos;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

void AllocateData(N_Vector v)
{
  N_VectorContent_Kokkos vc = NVEC_KOKKOS_CONTENT(v);

  vc->device_data = DeviceArrayView("device_data", NVEC_KOKKOS_MEMSIZE(v));
  vc->host_data = HostArrayView("host_data", NVEC_KOKKOS_MEMSIZE(v));

}

} // extern "C"
